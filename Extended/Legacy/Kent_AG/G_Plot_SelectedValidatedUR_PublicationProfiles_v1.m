%% ========================================================================
% F_Plot_SelectedValidatedUR_PublicationProfiles_v1_2.m
% ========================================================================
% Purpose
%   Publication-facing plotting for validated raw ultradian rhythm candidates.
%
%   This script does not regenerate wavelets or scalograms. It uses:
%     1) HSubSupported_PeriodMap.mat from Code C
%     2) WP_TS__*.mat handoff files from Code B
%     3) original raw behavioural Excel files
%
% Outputs
%   - Data-driven validated period cluster histograms
%   - DL/LD-only transition-aligned phase coherence curves
%   - 24 h phase coherence profiles
%   - 24 h ridge-power profiles
%   - 24 h period-targeted activity component profiles
%
% Conceptual notes
%   - Selective-HSub is used only as the validation gate upstream.
%   - The biological signal plotted here remains the Raw signal.
%   - The activity-component traces are bandpass-filtered around validated
%     period clusters and then z-scored per mouse/photoperiod.
%   - The y-axis is labelled as "Activity (z-scored)" for readability.
%
% MATLAB: written for R2025a/R2025b
% ========================================================================

clearvars; close all; clc;

%% ----------------------------- SETTINGS ---------------------------------
SCRIPT_NAME    = 'F_Plot_SelectedValidatedUR_PublicationProfiles_v1_2';
SCRIPT_VERSION = '1.2';
RUN_TIMESTAMP  = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

% Main bands for publication-focused plotting
BANDS_TO_USE = ["UR_1_3", "UR_3_6"];

% Cluster selection
CLUSTER_MODE = "auto";       % "auto" or "manual"

% Manual mode example:
% MANUAL_CLUSTERS = table(["UR_1_3";"UR_3_6"], [1.75;4.00], [2.25;4.50], ...
%     'VariableNames', {'BandName','PeriodLow_h','PeriodHigh_h'});
MANUAL_CLUSTERS = table(string.empty(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'BandName','PeriodLow_h','PeriodHigh_h'});

% Auto histogram peak detection
PERIOD_BIN_WIDTH_H          = 0.25;   % histogram bin width for validated period clusters
PERIOD_HIST_SMOOTH_SPAN     = 3;      % moving average span over histogram bins
PEAK_RELATIVE_HEIGHT_FRAC   = 0.35;   % cluster window expands until smoothed histogram falls below this fraction of peak
MIN_CANDIDATES_PER_CLUSTER  = 5;      % soft minimum. If all peaks fail, strongest available peak is retained.
MIN_MICE_PER_CLUSTER        = 3;      % soft minimum. If all peaks fail, strongest available peak is retained.
FILTER_PERIOD_PADDING_FRAC  = 0.10;   % 10 percent padding for activity-component bandpass filter

% Transition-aligned phase coherence
PERI_WINDOW_H     = 6.0;
REL_BIN_WIDTH_H   = 0.5;
MIN_PHASE_POINTS_PER_BIN = 10;

% Adaptive windows prevent L20/L22 windows from spanning too deeply into the neighbouring transition
ADAPTIVE_TRANSITION_WINDOWS = true;
TRANSITION_WINDOW_FRACTION  = 0.45;

% 24 h profiles
ZT_BIN_WIDTH_H = 0.25;
MIN_PHASE_POINTS_PER_ZT_BIN = 10;

% LL handling
LL_PHOTOPERIOD_VALUE = 24.0;
LL_PROJECTED_REFERENCE_PHOTOPERIOD_H = 22.0;   % previous LD condition before LL
LL_PROJECTED_DARK_SHADE_ALPHA = 0.10;

% Activity component filtering
FILTER_ORDER = 3;
REQUIRE_AT_LEAST_N_CYCLES_FOR_FILTER = 3;

% Outputs
SAVE_DPI = 600;
FIG_EXT = '.jpg';             % JPEG only. No TIFF output in this script.
FIG_FONT = 'Times New Roman';

% Plot appearance
INDIVIDUAL_LINE_COLOUR = [0.72 0.72 0.72];
MEAN_LINE_COLOUR       = [0.85 0.00 0.00];
TRANSITION_LINE_COLOUR = [0.05 0.05 0.05];
DARK_SHADE_COLOUR      = [0.82 0.82 0.82];
PROJECTED_SHADE_COLOUR = [0.72 0.72 0.72];

%% ----------------------------- FOLDERS ----------------------------------
fprintf('%s started at %s\n', SCRIPT_NAME, RUN_TIMESTAMP);

acrossFolder = uigetdir(pwd, 'Select AcrossPhotoperiod_Input folder containing WP_TS__*.mat');
if isequal(acrossFolder, 0), error('No AcrossPhotoperiod_Input folder selected.'); end

rawFolder = uigetdir(pwd, 'Select folder containing original raw behavioural Excel files');
if isequal(rawFolder, 0), error('No raw behavioural data folder selected.'); end

outputParent = uigetdir(pwd, 'Select output results folder for publication plots');
if isequal(outputParent, 0), error('No output results folder selected.'); end

outRoot = fullfile(outputParent, 'SelectedValidatedUR_PublicationProfiles');
dirs = struct();
dirs.Root       = outRoot;
dirs.Tables     = fullfile(outRoot, 'Tables');
dirs.FigRoot    = fullfile(outRoot, 'Figures');
dirs.ClusterFig = fullfile(dirs.FigRoot, '01_PeriodClusters');
dirs.PhaseRel   = fullfile(dirs.FigRoot, '02_PhaseCoherence_DL_LD');
dirs.PhaseZT    = fullfile(dirs.FigRoot, '03_24h_PhaseCoherence');
dirs.RidgePower = fullfile(dirs.FigRoot, '04_24h_RidgePower');
dirs.Activity   = fullfile(dirs.FigRoot, '05_24h_Activity_zscored');
dirs.Logs       = fullfile(outRoot, 'Logs');

ensure_dir(dirs.Root);
ensure_dir(dirs.Tables);
ensure_dir(dirs.FigRoot);
ensure_dir(dirs.ClusterFig);
ensure_dir(dirs.PhaseRel);
ensure_dir(dirs.PhaseZT);
ensure_dir(dirs.RidgePower);
ensure_dir(dirs.Activity);
ensure_dir(dirs.Logs);

logPath = fullfile(dirs.Logs, sprintf('%s_RunLog_%s.txt', SCRIPT_NAME, datestr(now,'yyyymmdd_HHMMSS')));
LOG = fopen(logPath, 'w');
cleanupObj = onCleanup(@() fclose_if_valid(LOG)); %#ok<NASGU>
log_line(LOG, '%s started at %s', SCRIPT_NAME, RUN_TIMESTAMP);
log_line(LOG, 'AcrossPhotoperiod_Input: %s', acrossFolder);
log_line(LOG, 'Raw behavioural data folder: %s', rawFolder);
log_line(LOG, 'Output root: %s', outRoot);

%% ----------------------------- LOAD VALIDATION MAP -----------------------
valPath = fullfile(acrossFolder, 'RawVsSelectiveHSub_PeriodValidation', 'HSubSupported_PeriodMap.mat');
if ~isfile(valPath)
    [vf, vp] = uigetfile('*.mat', 'Select HSubSupported_PeriodMap.mat', acrossFolder);
    if isequal(vf,0), error('No validation map selected.'); end
    valPath = fullfile(vp, vf);
end

CF = load_carryforward_table(valPath);
CF = standardise_carryforward_table(CF);
CF = CF(ismember(string(CF.BandName), BANDS_TO_USE), :);
CF = CF(string(CF.Phase) == "All", :);
CF = CF(isfinite(CF.RawPeriod_h) & CF.RawPeriod_h > 0, :);

log_line(LOG, 'Loaded carry-forward validated Raw candidates: %d', height(CF));
fprintf('Carry-forward validated Raw candidates used: %d\n', height(CF));

if height(CF) == 0
    error('No carry-forward validated Raw candidates were found for the selected bands and Phase == All.');
end

%% ----------------------------- LOAD RIDGE PHASE TABLES -------------------
tsFiles = dir(fullfile(acrossFolder, 'WP_TS__*.mat'));
if isempty(tsFiles)
    error('No WP_TS__*.mat files found in: %s', acrossFolder);
end

fprintf('Found %d WP_TS files.\n', numel(tsFiles));
log_line(LOG, 'Found %d WP_TS files.', numel(tsFiles));

RidgePhase = table();
for i = 1:numel(tsFiles)
    fpath = fullfile(tsFiles(i).folder, tsFiles(i).name);
    try
        S = load(fpath);
        T = extract_ridgephase_table(S, fpath);
        T = standardise_ridgephase_table(T);
        T = T(ismember(string(T.BandName), BANDS_TO_USE), :);
        T = T(string(T.Source) == "Raw", :);
        T = T(string(T.Phase) == "All", :);
        T = T(logical(T.ValidFlag), :);
        if isempty(RidgePhase) || width(RidgePhase) == 0
            RidgePhase = T;
        else
            RidgePhase = [RidgePhase; T]; %#ok<AGROW>
        end
        fprintf('Loaded %s (%d usable RidgePhase rows).\n', tsFiles(i).name, height(T));
        log_line(LOG, 'Loaded %s (%d usable RidgePhase rows).', tsFiles(i).name, height(T));
    catch ME
        warning('Skipping %s: %s', tsFiles(i).name, ME.message);
        log_line(LOG, 'WARNING: skipping %s: %s', tsFiles(i).name, ME.message);
    end
end

if height(RidgePhase) == 0
    error('No usable Raw Phase == All RidgePhase rows were loaded.');
end

% Keep only validated candidate IDs
validIDs = unique(string(CF.CandidateID));
RidgePhase = RidgePhase(ismember(string(RidgePhase.CandidateID), validIDs), :);
log_line(LOG, 'RidgePhase rows after validated candidate filter: %d', height(RidgePhase));
fprintf('RidgePhase rows after validated candidate filter: %d\n', height(RidgePhase));

if height(RidgePhase) == 0
    error('RidgePhase candidate IDs did not match the validation map. Check CandidateID/RawCandidateID consistency.');
end

%% ----------------------------- CLUSTER VALIDATED PERIODS -----------------
Settings = make_settings_table(SCRIPT_NAME, SCRIPT_VERSION, RUN_TIMESTAMP, acrossFolder, rawFolder, outRoot, ...
    BANDS_TO_USE, CLUSTER_MODE, PERIOD_BIN_WIDTH_H, FILTER_PERIOD_PADDING_FRAC, REL_BIN_WIDTH_H, ZT_BIN_WIDTH_H);

[ClusterSummary, ClusterMembership, AllCandidatePeriods] = build_period_clusters(CF, BANDS_TO_USE, CLUSTER_MODE, MANUAL_CLUSTERS, ...
    PERIOD_BIN_WIDTH_H, PERIOD_HIST_SMOOTH_SPAN, PEAK_RELATIVE_HEIGHT_FRAC, ...
    MIN_CANDIDATES_PER_CLUSTER, MIN_MICE_PER_CLUSTER, LOG);

if height(ClusterSummary) == 0
    error('No validated period clusters were detected. Consider switching CLUSTER_MODE to manual or lowering thresholds.');
end

fprintf('Detected/defined %d validated period clusters.\n', height(ClusterSummary));
log_line(LOG, 'Detected/defined %d validated period clusters.', height(ClusterSummary));

% Add cluster metadata to RidgePhase
RidgePhaseClustered = innerjoin(RidgePhase, ClusterMembership(:, {'CandidateID','ClusterID','ClusterRank','PeriodLow_h','PeriodHigh_h','PeriodCentre_h','FilterLow_h','FilterHigh_h'}), ...
    'Keys', 'CandidateID');

log_line(LOG, 'RidgePhase rows after cluster membership join: %d', height(RidgePhaseClustered));
fprintf('RidgePhase rows after cluster membership join: %d\n', height(RidgePhaseClustered));

%% ----------------------------- RAW FILE INDEX ----------------------------
RawFileIndex = index_raw_excel_files(rawFolder);
log_line(LOG, 'Indexed %d raw Excel files.', height(RawFileIndex));

%% ----------------------------- OUTPUT TABLES -----------------------------
xlsxOut = fullfile(dirs.Tables, 'SelectedValidatedUR_PublicationProfiles_Output.xlsx');
safe_delete_file(xlsxOut);

safe_writetable(Settings, xlsxOut, 'Settings');
safe_writetable(AllCandidatePeriods, xlsxOut, 'AllCandidatePeriods');
safe_writetable(ClusterSummary, xlsxOut, 'ClusterSummary');
safe_writetable(ClusterMembership, xlsxOut, 'ClusterMembership');
safe_writetable(RawFileIndex, xlsxOut, 'RawFileIndex');

PlotFiles = table(string.empty(0,1), string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    'VariableNames', {'PlotType','ClusterID','Photoperiod','FilePath'});

%% ----------------------------- PERIOD HISTOGRAMS -------------------------
fprintf('Plotting validated period cluster histograms...\n');
[PlotFiles_cluster] = plot_period_cluster_histograms(AllCandidatePeriods, ClusterSummary, dirs.ClusterFig, ...
    BANDS_TO_USE, PERIOD_BIN_WIDTH_H, SAVE_DPI, FIG_EXT, FIG_FONT);
PlotFiles = [PlotFiles; PlotFiles_cluster];

%% ----------------------------- PHASE COHERENCE ---------------------------
fprintf('Plotting DL/LD-only phase coherence curves...\n');
[PhaseRelTable, PlotFiles_rel] = plot_transition_phase_coherence(RidgePhaseClustered, ClusterSummary, dirs.PhaseRel, ...
    PERI_WINDOW_H, REL_BIN_WIDTH_H, MIN_PHASE_POINTS_PER_BIN, ...
    ADAPTIVE_TRANSITION_WINDOWS, TRANSITION_WINDOW_FRACTION, ...
    LL_PHOTOPERIOD_VALUE, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H, SAVE_DPI, FIG_EXT, FIG_FONT, ...
    MEAN_LINE_COLOUR, LOG);
PlotFiles = [PlotFiles; PlotFiles_rel];
safe_writetable(PhaseRelTable, xlsxOut, 'PhaseCoherence_DL_LD');

fprintf('Plotting 24 h phase coherence profiles...\n');
[PhaseZTTable, PlotFiles_phaseZT] = plot_zt_phase_coherence(RidgePhaseClustered, ClusterSummary, dirs.PhaseZT, ...
    ZT_BIN_WIDTH_H, MIN_PHASE_POINTS_PER_ZT_BIN, LL_PHOTOPERIOD_VALUE, ...
    LL_PROJECTED_REFERENCE_PHOTOPERIOD_H, SAVE_DPI, FIG_EXT, FIG_FONT, LOG);
PlotFiles = [PlotFiles; PlotFiles_phaseZT];
safe_writetable(PhaseZTTable, xlsxOut, 'PhaseCoherence_24h');

%% ----------------------------- RIDGE POWER PROFILES ----------------------
fprintf('Plotting 24 h ridge-power profiles...\n');
[RidgePowerZTTable, PlotFiles_rp] = plot_zt_ridge_power(RidgePhaseClustered, ClusterSummary, dirs.RidgePower, ...
    ZT_BIN_WIDTH_H, LL_PHOTOPERIOD_VALUE, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H, ...
    SAVE_DPI, FIG_EXT, FIG_FONT, INDIVIDUAL_LINE_COLOUR, MEAN_LINE_COLOUR, LOG);
PlotFiles = [PlotFiles; PlotFiles_rp];
safe_writetable(RidgePowerZTTable, xlsxOut, 'RidgePower_24h');

%% ----------------------------- ACTIVITY COMPONENT PROFILES ---------------
fprintf('Plotting 24 h period-targeted activity profiles...\n');
[ActivityZTTable, ActivityWarnings, PlotFiles_act] = plot_zt_activity_components(ClusterSummary, ClusterMembership, RawFileIndex, dirs.Activity, ...
    ZT_BIN_WIDTH_H, FILTER_PERIOD_PADDING_FRAC, FILTER_ORDER, REQUIRE_AT_LEAST_N_CYCLES_FOR_FILTER, ...
    LL_PHOTOPERIOD_VALUE, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H, SAVE_DPI, FIG_EXT, FIG_FONT, ...
    INDIVIDUAL_LINE_COLOUR, MEAN_LINE_COLOUR, LOG);
PlotFiles = [PlotFiles; PlotFiles_act];
safe_writetable(ActivityZTTable, xlsxOut, 'ActivityComponent_24h');
safe_writetable(ActivityWarnings, xlsxOut, 'ActivityWarnings');

safe_writetable(PlotFiles, xlsxOut, 'PlotFiles');

matOut = fullfile(dirs.Tables, 'SelectedValidatedUR_PublicationProfiles_Output.mat');
save(matOut, 'Settings', 'AllCandidatePeriods', 'ClusterSummary', 'ClusterMembership', ...
    'PhaseRelTable', 'PhaseZTTable', 'RidgePowerZTTable', 'ActivityZTTable', ...
    'ActivityWarnings', 'PlotFiles', '-v7.3');

fprintf('\nComplete.\n');
fprintf('Tables written to:\n  %s\n', xlsxOut);
fprintf('Figures written under:\n  %s\n', dirs.FigRoot);
log_line(LOG, 'Complete. Tables: %s', xlsxOut);
log_line(LOG, 'Figures: %s', dirs.FigRoot);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function CF = load_carryforward_table(valPath)
    S = load(valPath);
    CF = table();

    if isfield(S, 'validationMap')
        vm = S.validationMap;
        if isstruct(vm) && isfield(vm, 'CarryForward_Periods')
            CF = vm.CarryForward_Periods;
        elseif istable(vm)
            CF = vm;
        end
    end

    if isempty(CF)
        fns = fieldnames(S);
        for i = 1:numel(fns)
            obj = S.(fns{i});
            if istable(obj) && any(strcmp(obj.Properties.VariableNames, 'RawCandidateID'))
                CF = obj;
                break;
            end
        end
    end

    if isempty(CF) || ~istable(CF)
        error('Could not find CarryForward_Periods table in validation map: %s', valPath);
    end
end

function CF = standardise_carryforward_table(CF)
    required = {'File','SignalID','Photoperiod_h','Phase','BandName','RawCandidateID','RawPeriod_h'};
    for i = 1:numel(required)
        if ~ismember(required{i}, CF.Properties.VariableNames)
            error('CarryForward table is missing required column: %s', required{i});
        end
    end

    if ismember('CarryForward', CF.Properties.VariableNames)
        CF = CF(to_logical(CF.CarryForward), :);
    end

    CF.File = string(CF.File);
    CF.SignalID = string(CF.SignalID);
    if ismember('ConditionParsed', CF.Properties.VariableNames)
        CF.ConditionParsed = string(CF.ConditionParsed);
    else
        CF.ConditionParsed = strings(height(CF),1);
    end
    CF.Photoperiod_h = double(CF.Photoperiod_h);
    CF.Phase = string(CF.Phase);
    CF.BandName = string(CF.BandName);
    CF.RawCandidateID = string(CF.RawCandidateID);
    CF.CandidateID = string(CF.RawCandidateID);
    CF.RawPeriod_h = double(CF.RawPeriod_h);

    optionalNum = {'RawIQR_h','PrimaryHSubPeriod_h','PeriodDiffPercent','RawRidgeCoverageFrac','RawCOIValidFrac'};
    for i = 1:numel(optionalNum)
        v = optionalNum{i};
        if ismember(v, CF.Properties.VariableNames)
            CF.(v) = double(CF.(v));
        else
            CF.(v) = NaN(height(CF),1);
        end
    end

    optionalStr = {'FullLadderSensitivityStatus','FinalValidationClass'};
    for i = 1:numel(optionalStr)
        v = optionalStr{i};
        if ismember(v, CF.Properties.VariableNames)
            CF.(v) = string(CF.(v));
        else
            CF.(v) = strings(height(CF),1);
        end
    end
end

function T = extract_ridgephase_table(S, sourceFile)
    T = table();
    if isfield(S, 'pkgTS') && isstruct(S.pkgTS)
        pkgTS = S.pkgTS;
        if isfield(pkgTS, 'tables') && isstruct(pkgTS.tables) && isfield(pkgTS.tables, 'RidgePhase_Long')
            T = pkgTS.tables.RidgePhase_Long;
        end
    end

    if isempty(T)
        fns = fieldnames(S);
        for i = 1:numel(fns)
            obj = S.(fns{i});
            if istable(obj) && all(ismember({'CandidateID','RidgePhase_rad','RidgePower_log10'}, obj.Properties.VariableNames))
                T = obj;
                break;
            end
        end
    end

    if isempty(T) || ~istable(T)
        error('No RidgePhase_Long table found in %s', sourceFile);
    end
end

function T = standardise_ridgephase_table(T)
    required = {'File','SignalID','Source','Photoperiod_h','BandName','CandidateID', ...
                'Time_days','ZT_hours','Phase','RidgePeriod_h','RidgePower_log10','RidgePhase_rad','ValidFlag'};
    for i = 1:numel(required)
        if ~ismember(required{i}, T.Properties.VariableNames)
            error('RidgePhase_Long is missing required column: %s', required{i});
        end
    end

    strVars = {'File','SignalID','Source','BandName','CandidateID','Phase'};
    for i = 1:numel(strVars)
        T.(strVars{i}) = string(T.(strVars{i}));
    end
    if ismember('ConditionParsed', T.Properties.VariableNames)
        T.ConditionParsed = string(T.ConditionParsed);
    else
        T.ConditionParsed = strings(height(T),1);
    end
    if ismember('LightStateValue', T.Properties.VariableNames)
        T.LightStateValue = string(T.LightStateValue);
    else
        T.LightStateValue = strings(height(T),1);
    end
    if ismember('HSubResidualMode', T.Properties.VariableNames)
        T.HSubResidualMode = string(T.HSubResidualMode);
    else
        T.HSubResidualMode = strings(height(T),1);
    end

    numVars = {'Photoperiod_h','Time_days','ZT_hours','RidgePeriod_h','RidgePower_log10','RidgePhase_rad'};
    for i = 1:numel(numVars)
        T.(numVars{i}) = double(T.(numVars{i}));
    end
    T.ValidFlag = to_logical(T.ValidFlag);
end

function Settings = make_settings_table(scriptName, scriptVersion, runTimestamp, acrossFolder, rawFolder, outRoot, bands, clusterMode, binWidth, padFrac, relBin, ztBin)
    names = ["ScriptName"; "ScriptVersion"; "RunTimestamp"; "AcrossPhotoperiodInput"; "RawDataFolder"; "OutputRoot"; ...
        "Bands"; "ClusterMode"; "PeriodHistogramBinWidth_h"; "ActivityFilterPaddingFrac"; ...
        "RelativePhaseCoherenceBinWidth_h"; "ZTProfileBinWidth_h"; "ActivityYAxisLabel"; "TIFFOutput"];
    vals = [string(scriptName); string(scriptVersion); string(runTimestamp); string(acrossFolder); string(rawFolder); string(outRoot); ...
        strjoin(bands, ", "); string(clusterMode); string(binWidth); string(padFrac); ...
        string(relBin); string(ztBin); "Activity (z-scored)"; "Disabled"];
    Settings = table(names, vals, 'VariableNames', {'Setting','Value'});
end

function [ClusterSummary, ClusterMembership, AllCand] = build_period_clusters(CF, bands, mode, manualClusters, binWidth, smoothSpan, relHeightFrac, minCand, minMice, LOG)
    AllCand = table();
    AllCand.CandidateID = string(CF.CandidateID);
    AllCand.File = string(CF.File);
    AllCand.SignalID = string(CF.SignalID);
    AllCand.ConditionParsed = string(CF.ConditionParsed);
    AllCand.Photoperiod_h = double(CF.Photoperiod_h);
    AllCand.BandName = string(CF.BandName);
    AllCand.RawPeriod_h = double(CF.RawPeriod_h);
    AllCand.RawIQR_h = double(CF.RawIQR_h);
    AllCand.PrimaryHSubPeriod_h = double(CF.PrimaryHSubPeriod_h);
    AllCand.PeriodDiffPercent = double(CF.PeriodDiffPercent);

    csRows = {};
    cmRows = {};
    csHdr = {'BandName','ClusterID','ClusterRank','PeriodLow_h','PeriodHigh_h','PeriodCentre_h', ...
        'FilterLow_h','FilterHigh_h','CandidateCount','MouseCount','PhotoperiodCount','PhotoperiodsRepresented', ...
        'MedianRawPeriod_h','IQRRawPeriod_h','ClusterMode','SelectionNote'};
    cmHdr = {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','RawPeriod_h', ...
        'ClusterID','ClusterRank','PeriodLow_h','PeriodHigh_h','PeriodCentre_h','FilterLow_h','FilterHigh_h'};

    for b = 1:numel(bands)
        band = string(bands(b));
        B = AllCand(string(AllCand.BandName) == band, :);
        if height(B) == 0
            log_line(LOG, 'No candidates for band %s.', band);
            continue;
        end

        windows = [];
        notes = strings(0,1);

        if string(mode) == "manual"
            M = manualClusters(string(manualClusters.BandName) == band, :);
            for i = 1:height(M)
                windows(end+1,:) = [double(M.PeriodLow_h(i)), double(M.PeriodHigh_h(i))]; %#ok<AGROW>
                notes(end+1,1) = "Manual user-defined window"; %#ok<AGROW>
            end
        else
            [windows, notes] = detect_histogram_peak_windows(B.RawPeriod_h, band, binWidth, smoothSpan, relHeightFrac, minCand, minMice, B.SignalID);
        end

        if isempty(windows)
            log_line(LOG, 'No cluster windows retained for %s.', band);
            continue;
        end

        % Build the ranking table using numeric vectors rather than growing an
        % empty table by dot-indexing. Some MATLAB releases initialise empty
        % table variables as cell-like during this pattern and then fail when
        % numeric values are assigned.
        nWin = size(windows,1);
        metric_Row = NaN(nWin,1);
        metric_CandidateCount = NaN(nWin,1);
        metric_MouseCount = NaN(nWin,1);
        metric_PhotoperiodCount = NaN(nWin,1);
        metric_MedianRawPeriod_h = NaN(nWin,1);

        for w = 1:nWin
            low = windows(w,1);
            high = windows(w,2);
            idx = B.RawPeriod_h >= low & B.RawPeriod_h <= high;
            Bi = B(idx,:);

            metric_Row(w) = w;
            metric_CandidateCount(w) = height(Bi);
            metric_MouseCount(w) = count_unique(Bi.SignalID);
            metric_PhotoperiodCount(w) = count_unique(Bi.Photoperiod_h);
            metric_MedianRawPeriod_h(w) = median(Bi.RawPeriod_h, 'omitnan');
        end

        metrics = table(metric_Row, metric_CandidateCount, metric_MouseCount, ...
            metric_PhotoperiodCount, metric_MedianRawPeriod_h, ...
            'VariableNames', {'WindowIndex','CandidateCount','MouseCount','PhotoperiodCount','MedianRawPeriod_h'});

        metrics = sortrows(metrics, {'CandidateCount','MouseCount','PhotoperiodCount','MedianRawPeriod_h'}, {'descend','descend','descend','ascend'});

        for rank = 1:height(metrics)
            w = metrics.WindowIndex(rank);
            low = windows(w,1);
            high = windows(w,2);
            centre = mean([low high], 'omitnan');
            [fLow, fHigh] = padded_period_window(low, high, 0.10);
            idx = B.RawPeriod_h >= low & B.RawPeriod_h <= high;
            Bi = B(idx,:);

            if isempty(Bi), continue; end

            clusterID = sprintf('%s_C%02d_%s_%sh', sanitise_filename(band), rank, safe_period(low), safe_period(high));
            photoVals = unique(Bi.Photoperiod_h(isfinite(Bi.Photoperiod_h)));
            photoStr = strjoin(arrayfun(@(x) photo_label(x), photoVals, 'UniformOutput', false), ', ');

            csRows(end+1,:) = {band, string(clusterID), rank, low, high, centre, fLow, fHigh, ...
                height(Bi), count_unique(Bi.SignalID), count_unique(Bi.Photoperiod_h), string(photoStr), ...
                median(Bi.RawPeriod_h,'omitnan'), local_iqr(Bi.RawPeriod_h), string(mode), string(notes(w))}; %#ok<AGROW>

            for r = 1:height(Bi)
                cmRows(end+1,:) = {string(Bi.CandidateID(r)), string(Bi.File(r)), string(Bi.SignalID(r)), string(Bi.ConditionParsed(r)), ...
                    double(Bi.Photoperiod_h(r)), string(Bi.BandName(r)), double(Bi.RawPeriod_h(r)), ...
                    string(clusterID), rank, low, high, centre, fLow, fHigh}; %#ok<AGROW>
            end
        end
    end

    ClusterSummary = rows_to_table(csRows, csHdr);
    ClusterMembership = rows_to_table(cmRows, cmHdr);
end

function [windows, notes] = detect_histogram_peak_windows(periods, band, binWidth, smoothSpan, relHeightFrac, minCand, minMice, signalIDs)
    periods = periods(:);
    periods = periods(isfinite(periods) & periods > 0);
    signalIDs = string(signalIDs(:));

    [bandLow, bandHigh] = parse_band_bounds_h(band);
    if ~isfinite(bandLow), bandLow = floor(min(periods)/binWidth)*binWidth; end
    if ~isfinite(bandHigh), bandHigh = ceil(max(periods)/binWidth)*binWidth; end

    edges = bandLow:binWidth:bandHigh;
    if numel(edges) < 3
        edges = floor(min(periods)/binWidth)*binWidth:binWidth:ceil(max(periods)/binWidth)*binWidth;
    end
    if numel(edges) < 3
        windows = [];
        notes = strings(0,1);
        return;
    end

    counts = histcounts(periods, edges);
    if smoothSpan > 1
        sm = movmean(counts, smoothSpan);
    else
        sm = counts;
    end

    peakIdx = local_peak_indices(sm);
    if isempty(peakIdx)
        [~, peakIdx] = max(sm);
    end

    provisional = [];
    pnotes = strings(0,1);
    for p = 1:numel(peakIdx)
        pk = peakIdx(p);
        if sm(pk) <= 0, continue; end
        thresh = max(1, relHeightFrac * sm(pk));

        L = pk;
        while L > 1 && sm(L-1) >= thresh && counts(L-1) > 0
            L = L - 1;
        end

        R = pk;
        while R < numel(sm) && sm(R+1) >= thresh && counts(R+1) > 0
            R = R + 1;
        end

        low = edges(L);
        high = edges(R+1);
        idx = periods >= low & periods <= high;
        nCand = sum(idx);
        nMice = count_unique(signalIDs(idx));
        note = sprintf('Auto histogram peak at %.3g h; candidates=%d; mice=%d', edges(pk)+binWidth/2, nCand, nMice);

        provisional(end+1,:) = [low high nCand nMice sm(pk)]; %#ok<AGROW>
        pnotes(end+1,1) = string(note); %#ok<AGROW>
    end

    if isempty(provisional)
        windows = [];
        notes = strings(0,1);
        return;
    end

    keep = provisional(:,3) >= minCand & provisional(:,4) >= minMice;
    if ~any(keep)
        [~, best] = max(provisional(:,5));
        keep(best) = true;
        pnotes(best) = pnotes(best) + "; retained as strongest available peak despite soft thresholds";
    end

    provisional = provisional(keep,:);
    pnotes = pnotes(keep);

    % Merge overlapping windows within a band
    [~, ord] = sort(provisional(:,1));
    provisional = provisional(ord,:);
    pnotes = pnotes(ord);

    merged = [];
    mergedNotes = strings(0,1);
    for i = 1:size(provisional,1)
        if isempty(merged)
            merged = provisional(i,:);
            mergedNotes(1,1) = pnotes(i);
        else
            if provisional(i,1) <= merged(end,2)
                merged(end,2) = max(merged(end,2), provisional(i,2));
                merged(end,3) = NaN;
                merged(end,4) = NaN;
                merged(end,5) = max(merged(end,5), provisional(i,5));
                mergedNotes(end) = mergedNotes(end) + "; merged with adjacent/overlapping peak";
            else
                merged(end+1,:) = provisional(i,:); %#ok<AGROW>
                mergedNotes(end+1,1) = pnotes(i); %#ok<AGROW>
            end
        end
    end

    windows = merged(:,1:2);
    notes = mergedNotes;
end

function idx = local_peak_indices(y)
    y = double(y(:));
    idx = [];
    n = numel(y);
    if n == 0, return; end
    for i = 1:n
        leftOK = i == 1 || y(i) >= y(i-1);
        rightOK = i == n || y(i) >= y(i+1);
        strict = (i == 1 || y(i) > y(i-1)) || (i == n || y(i) > y(i+1));
        if leftOK && rightOK && strict && y(i) > 0
            idx(end+1,1) = i; %#ok<AGROW>
        end
    end
end

function [lowPad, highPad] = padded_period_window(low, high, padFrac)
    lowPad = low * (1 - padFrac);
    highPad = high * (1 + padFrac);
    lowPad = max(lowPad, eps);
    highPad = max(highPad, lowPad + eps);
end

function [PlotFiles] = plot_period_cluster_histograms(AllCand, ClusterSummary, figDir, bands, binWidth, dpi, ext, figFont)
    PlotFiles = empty_plotfiles();
    for b = 1:numel(bands)
        band = string(bands(b));
        B = AllCand(string(AllCand.BandName) == band, :);
        C = ClusterSummary(string(ClusterSummary.BandName) == band, :);
        if height(B) == 0, continue; end

        out = fullfile(figDir, sprintf('PeriodCluster_Histogram_Pooled_%s%s', sanitise_filename(band), ext));
        make_one_period_histogram(B.RawPeriod_h, C, sprintf('Validated periods | pooled | %s', band), out, binWidth, dpi, figFont);
        PlotFiles = add_plotfile(PlotFiles, "PeriodClusterHistogram_Pooled", "", "Pooled", out);

        photos = unique(B.Photoperiod_h(isfinite(B.Photoperiod_h)));
        for p = 1:numel(photos)
            photo = photos(p);
            Bp = B(B.Photoperiod_h == photo, :);
            out = fullfile(figDir, sprintf('PeriodCluster_Histogram_%s_%s%s', photo_label(photo), sanitise_filename(band), ext));
            make_one_period_histogram(Bp.RawPeriod_h, C, sprintf('Validated periods | %s | %s', photo_label(photo), band), out, binWidth, dpi, figFont);
            PlotFiles = add_plotfile(PlotFiles, "PeriodClusterHistogram_ByPhotoperiod", "", photo_label(photo), out);
        end
    end
end

function make_one_period_histogram(periods, C, ttl, out, binWidth, dpi, figFont)
    periods = periods(isfinite(periods));
    if isempty(periods), return; end
    lo = floor(min(periods)/binWidth)*binWidth;
    hi = ceil(max(periods)/binWidth)*binWidth;
    if hi <= lo, hi = lo + binWidth; end
    edges = lo:binWidth:hi;
    if numel(edges) < 3, edges = [lo lo+binWidth lo+2*binWidth]; end

    f = figure('Color','w','Position',[100 100 900 520]);
    ax = axes(f); hold(ax,'on');
    histogram(ax, periods, edges, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'k');
    ylabel(ax, 'Validated candidates', 'FontWeight','bold');
    xlabel(ax, 'Median validated raw ridge period (h)', 'FontWeight','bold');
    title(ax, ttl, 'Interpreter','none', 'FontWeight','bold');

    yl = ylim(ax);
    for i = 1:height(C)
        x1 = C.PeriodLow_h(i);
        x2 = C.PeriodHigh_h(i);
        ph = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [1.0 0.75 0.75], ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
        try, uistack(ph, 'bottom'); catch, end
        xline(ax, C.PeriodCentre_h(i), ':', sprintf('C%d', C.ClusterRank(i)), ...
            'LineWidth', 1.2, 'LabelOrientation','horizontal', 'HandleVisibility','off');
    end

    ylim(ax, yl);
    box(ax,'off');
    set(ax,'TickDir','out','FontName',figFont,'FontSize',14,'LineWidth',1.1);
    exportgraphics(f, out, 'Resolution', dpi);
    close(f);
end

function [PhaseRelTable, PlotFiles] = plot_transition_phase_coherence(RP, ClusterSummary, figDir, periWindowH, binWidth, minN, adaptive, winFrac, llValue, llRefPhoto, dpi, ext, figFont, meanColour, LOG)
    PlotFiles = empty_plotfiles();
    rows = {};
    hdr = {'ClusterID','ClusterRank','BandName','Photoperiod_h','TransitionType','RelBinCenter_h','RelBinStart_h','RelBinEnd_h', ...
        'N_PhaseObs','N_Mice','N_Candidates','R','MeanPhase_rad','CircularSD_rad','MeanRidgePeriod_h','MeanRidgePower_log10','PassN'};

    combos = unique(RP(:, {'ClusterID','Photoperiod_h'}), 'rows', 'stable');
    for i = 1:height(combos)
        clusterID = string(combos.ClusterID(i));
        photo = double(combos.Photoperiod_h(i));
        C = ClusterSummary(string(ClusterSummary.ClusterID) == clusterID, :);
        if isempty(C), continue; end

        P = RP(string(RP.ClusterID) == clusterID & RP.Photoperiod_h == photo, :);
        if height(P) == 0, continue; end

        isLL = photo >= llValue || photo >= 24;
        if isLL
            transTypes = ["ProjectedDL_LL", "ProjectedLD_LL"];
            baseTypes = ["DL", "LD"];
            zts = [0, llRefPhoto];
            photoForWindows = llRefPhoto;
        else
            transTypes = ["DL", "LD"];
            baseTypes = ["DL", "LD"];
            zts = [0, photo];
            photoForWindows = photo;
        end

        f = figure('Color','w','Position',[100 100 950 520]);
        ax = axes(f); hold(ax,'on');
        lineColours = [0.00 0.45 0.74; 0.85 0.33 0.10];

        for t = 1:numel(transTypes)
            [preWin, postWin] = transition_window_limits(photoForWindows, baseTypes(t), periWindowH, adaptive, winFrac);
            if preWin <= 0 || postWin <= 0, continue; end
            TL = collect_relative_phase_rows(P, transTypes(t), zts(t), preWin, postWin);
            if isempty(TL) || height(TL) == 0, continue; end

            edges = -preWin:binWidth:postWin;
            if edges(end) < postWin, edges(end+1) = postWin; end
            centers = edges(1:end-1) + diff(edges)/2;
            Rvals = NaN(numel(centers),1);

            for b = 1:numel(centers)
                idx = TL.RelativeTime_h >= edges(b) & TL.RelativeTime_h < edges(b+1);
                if b == numel(centers)
                    idx = TL.RelativeTime_h >= edges(b) & TL.RelativeTime_h <= edges(b+1);
                end
                theta = TL.RidgePhase_rad(idx);
                theta = theta(isfinite(theta));
                [R, mu, csd] = circ_summary(theta);
                n = numel(theta);
                pass = n >= minN;
                Rvals(b) = R;
                rows(end+1,:) = {clusterID, C.ClusterRank(1), string(C.BandName(1)), photo, transTypes(t), centers(b), edges(b), edges(b+1), ...
                    n, count_unique(TL.SignalID(idx)), count_unique(TL.CandidateID(idx)), R, mu, csd, ...
                    mean(TL.RidgePeriod_h(idx),'omitnan'), mean(TL.RidgePower_log10(idx),'omitnan'), pass}; %#ok<AGROW>
            end

            plot(ax, centers, Rvals, '-o', 'LineWidth', 1.8, 'MarkerSize', 4, ...
                'Color', lineColours(t,:), 'MarkerFaceColor', 'w', 'DisplayName', transition_display_name(transTypes(t)));
        end

        xline(ax, 0, ':', 'Transition', 'LineWidth', 1.2, 'Color', [0.1 0.1 0.1], 'HandleVisibility','off');
        xlabel(ax, 'Time relative to transition (h)', 'FontWeight','bold');
        ylabel(ax, 'Phase coherence, R', 'FontWeight','bold');
        title(ax, sprintf('DL/LD phase coherence | %s | %s | C%d %.2g-%.2g h', ...
            photo_label(photo), string(C.BandName(1)), C.ClusterRank(1), C.PeriodLow_h(1), C.PeriodHigh_h(1)), ...
            'Interpreter','none', 'FontWeight','bold');
        ylim(ax, [0 1]);
        box(ax,'off');
        set(ax,'TickDir','out','FontName',figFont,'FontSize',14,'LineWidth',1.1);
        legend(ax, 'Location','eastoutside', 'Interpreter','none');

        out = fullfile(figDir, sprintf('PhaseCoherence_DL_LD_%s_%s%s', photo_label(photo), sanitise_filename(clusterID), ext));
        exportgraphics(f, out, 'Resolution', dpi);
        close(f);
        PlotFiles = add_plotfile(PlotFiles, "PhaseCoherence_DL_LD", clusterID, photo_label(photo), out);
        log_line(LOG, 'Wrote %s', out);
    end

    PhaseRelTable = rows_to_table(rows, hdr);
end

function TL = collect_relative_phase_rows(P, transitionType, transitionZT, preWin, postWin)
    rows = {};
    hdr = {'File','SignalID','CandidateID','Photoperiod_h','BandName','ClusterID','TransitionType','TransitionZT_h','TransitionDay', ...
        'RelativeTime_h','Time_days','ZT_hours','RidgePhase_rad','RidgePeriod_h','RidgePower_log10'};

    candIDs = unique(string(P.CandidateID), 'stable');
    for c = 1:numel(candIDs)
        Pc = P(string(P.CandidateID) == candIDs(c), :);
        if isempty(Pc), continue; end
        tH = double(Pc.Time_days) * 24;
        if all(~isfinite(tH)), continue; end
        tMin = min(tH, [], 'omitnan');
        tMax = max(tH, [], 'omitnan');
        dayStart = floor((tMin - postWin) / 24) - 1;
        dayEnd   = ceil((tMax + preWin) / 24) + 1;

        for d = dayStart:dayEnd
            evAbs = d*24 + transitionZT;
            rel = tH - evAbs;
            idx = rel >= -preWin & rel <= postWin & isfinite(Pc.RidgePhase_rad);
            if ~any(idx), continue; end
            Pi = Pc(idx,:);
            ri = rel(idx);
            for r = 1:height(Pi)
                rows(end+1,:) = {string(Pi.File(r)), string(Pi.SignalID(r)), string(Pi.CandidateID(r)), double(Pi.Photoperiod_h(r)), ...
                    string(Pi.BandName(r)), string(Pi.ClusterID(r)), string(transitionType), double(transitionZT), d, ...
                    ri(r), double(Pi.Time_days(r)), double(Pi.ZT_hours(r)), double(Pi.RidgePhase_rad(r)), ...
                    double(Pi.RidgePeriod_h(r)), double(Pi.RidgePower_log10(r))}; %#ok<AGROW>
            end
        end
    end
    TL = rows_to_table(rows, hdr);
end

function [PhaseZTTable, PlotFiles] = plot_zt_phase_coherence(RP, ClusterSummary, figDir, ztBinWidth, minN, llValue, llRefPhoto, dpi, ext, figFont, LOG)
    PlotFiles = empty_plotfiles();
    rows = {};
    hdr = {'ClusterID','ClusterRank','BandName','Photoperiod_h','ZTBinCenter_h','ZTBinStart_h','ZTBinEnd_h', ...
        'N_PhaseObs','N_Mice','N_Candidates','R','MeanPhase_rad','CircularSD_rad','MeanRidgePeriod_h','MeanRidgePower_log10','PassN'};

    edges = 0:ztBinWidth:24;
    centers = edges(1:end-1) + diff(edges)/2;

    combos = unique(RP(:, {'ClusterID','Photoperiod_h'}), 'rows', 'stable');
    for i = 1:height(combos)
        clusterID = string(combos.ClusterID(i));
        photo = double(combos.Photoperiod_h(i));
        C = ClusterSummary(string(ClusterSummary.ClusterID) == clusterID, :);
        if isempty(C), continue; end
        P = RP(string(RP.ClusterID) == clusterID & RP.Photoperiod_h == photo, :);
        if height(P) == 0, continue; end

        Rvals = NaN(numel(centers),1);
        for b = 1:numel(centers)
            idx = P.ZT_hours >= edges(b) & P.ZT_hours < edges(b+1);
            if b == numel(centers)
                idx = P.ZT_hours >= edges(b) & P.ZT_hours <= edges(b+1);
            end
            theta = P.RidgePhase_rad(idx);
            theta = theta(isfinite(theta));
            [R, mu, csd] = circ_summary(theta);
            n = numel(theta);
            pass = n >= minN;
            Rvals(b) = R;
            rows(end+1,:) = {clusterID, C.ClusterRank(1), string(C.BandName(1)), photo, centers(b), edges(b), edges(b+1), ...
                n, count_unique(P.SignalID(idx)), count_unique(P.CandidateID(idx)), R, mu, csd, ...
                mean(P.RidgePeriod_h(idx),'omitnan'), mean(P.RidgePower_log10(idx),'omitnan'), pass}; %#ok<AGROW>
        end

        f = figure('Color','w','Position',[100 100 950 520]);
        ax = axes(f); hold(ax,'on');
        plot(ax, centers, Rvals, '-o', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', [0.00 0.45 0.74], 'MarkerFaceColor','w');
        add_light_dark_annotations(ax, photo, llValue, llRefPhoto, [0.82 0.82 0.82], [0.72 0.72 0.72], [0.05 0.05 0.05], 0.10);
        xlabel(ax, 'ZT / projected ZT (h)', 'FontWeight','bold');
        ylabel(ax, 'Phase coherence, R', 'FontWeight','bold');
        title(ax, sprintf('24 h phase coherence | %s | %s | C%d %.2g-%.2g h', ...
            photo_label(photo), string(C.BandName(1)), C.ClusterRank(1), C.PeriodLow_h(1), C.PeriodHigh_h(1)), ...
            'Interpreter','none', 'FontWeight','bold');
        xlim(ax, [0 24]); ylim(ax, [0 1]);
        box(ax,'off');
        set(ax,'TickDir','out','FontName',figFont,'FontSize',14,'LineWidth',1.1);

        out = fullfile(figDir, sprintf('PhaseCoherence_24h_%s_%s%s', photo_label(photo), sanitise_filename(clusterID), ext));
        exportgraphics(f, out, 'Resolution', dpi);
        close(f);
        PlotFiles = add_plotfile(PlotFiles, "PhaseCoherence_24h", clusterID, photo_label(photo), out);
        log_line(LOG, 'Wrote %s', out);
    end

    PhaseZTTable = rows_to_table(rows, hdr);
end

function [RidgePowerZTTable, PlotFiles] = plot_zt_ridge_power(RP, ClusterSummary, figDir, ztBinWidth, llValue, llRefPhoto, dpi, ext, figFont, indivColour, meanColour, LOG)
    PlotFiles = empty_plotfiles();
    rows = {};
    hdr = {'ClusterID','ClusterRank','BandName','Photoperiod_h','SignalID','CandidateID','ZTBinCenter_h','ZTBinStart_h','ZTBinEnd_h', ...
        'MeanRidgePower_log10','N_Obs'};

    edges = 0:ztBinWidth:24;
    centers = edges(1:end-1) + diff(edges)/2;

    combos = unique(RP(:, {'ClusterID','Photoperiod_h'}), 'rows', 'stable');
    for i = 1:height(combos)
        clusterID = string(combos.ClusterID(i));
        photo = double(combos.Photoperiod_h(i));
        C = ClusterSummary(string(ClusterSummary.ClusterID) == clusterID, :);
        if isempty(C), continue; end
        P = RP(string(RP.ClusterID) == clusterID & RP.Photoperiod_h == photo, :);
        if height(P) == 0, continue; end

        candIDs = unique(string(P.CandidateID), 'stable');
        Y = NaN(numel(candIDs), numel(centers));

        for c = 1:numel(candIDs)
            Pc = P(string(P.CandidateID) == candIDs(c), :);
            sig = robust_first_string(Pc.SignalID);
            for b = 1:numel(centers)
                idx = Pc.ZT_hours >= edges(b) & Pc.ZT_hours < edges(b+1);
                if b == numel(centers)
                    idx = Pc.ZT_hours >= edges(b) & Pc.ZT_hours <= edges(b+1);
                end
                val = mean(Pc.RidgePower_log10(idx), 'omitnan');
                Y(c,b) = val;
                rows(end+1,:) = {clusterID, C.ClusterRank(1), string(C.BandName(1)), photo, sig, candIDs(c), centers(b), edges(b), edges(b+1), val, sum(idx)}; %#ok<AGROW>
            end
        end

        f = figure('Color','w','Position',[100 100 950 520]);
        ax = axes(f); hold(ax,'on');
        for c = 1:size(Y,1)
            plot_background_line(ax, centers, Y(c,:), indivColour);
        end
        meanY = mean(Y, 1, 'omitnan');
        plot(ax, centers, meanY, '-', 'LineWidth', 2.5, 'Color', meanColour);
        add_light_dark_annotations(ax, photo, llValue, llRefPhoto, [0.82 0.82 0.82], [0.72 0.72 0.72], [0.05 0.05 0.05], 0.10);
        xlabel(ax, 'ZT / projected ZT (h)', 'FontWeight','bold');
        ylabel(ax, 'Ridge power (log_{10})', 'FontWeight','bold');
        title(ax, sprintf('24 h ridge power | %s | %s | C%d %.2g-%.2g h', ...
            photo_label(photo), string(C.BandName(1)), C.ClusterRank(1), C.PeriodLow_h(1), C.PeriodHigh_h(1)), ...
            'Interpreter','none', 'FontWeight','bold');
        xlim(ax, [0 24]);
        box(ax,'off');
        set(ax,'TickDir','out','FontName',figFont,'FontSize',14,'LineWidth',1.1);

        out = fullfile(figDir, sprintf('RidgePower_24h_%s_%s%s', photo_label(photo), sanitise_filename(clusterID), ext));
        exportgraphics(f, out, 'Resolution', dpi);
        close(f);
        PlotFiles = add_plotfile(PlotFiles, "RidgePower_24h", clusterID, photo_label(photo), out);
        log_line(LOG, 'Wrote %s', out);
    end

    RidgePowerZTTable = rows_to_table(rows, hdr);
end

function [ActivityZTTable, Warnings, PlotFiles] = plot_zt_activity_components(ClusterSummary, ClusterMembership, RawFileIndex, figDir, ztBinWidth, padFrac, filterOrder, minCycles, llValue, llRefPhoto, dpi, ext, figFont, indivColour, meanColour, LOG)
    PlotFiles = empty_plotfiles();
    rows = {};
    warnRows = {};
    hdr = {'ClusterID','ClusterRank','BandName','Photoperiod_h','File','SignalID','ZTBinCenter_h','ZTBinStart_h','ZTBinEnd_h', ...
        'Activity_zscored','N_Obs','FilterLow_h','FilterHigh_h'};
    whdr = {'ClusterID','Photoperiod_h','File','SignalID','Warning'};

    edges = 0:ztBinWidth:24;
    centers = edges(1:end-1) + diff(edges)/2;

    for ci = 1:height(ClusterSummary)
        C = ClusterSummary(ci,:);
        clusterID = string(C.ClusterID);
        photos = unique(ClusterMembership.Photoperiod_h(string(ClusterMembership.ClusterID)==clusterID));
        photos = photos(isfinite(photos));

        for p = 1:numel(photos)
            photo = photos(p);
            M = ClusterMembership(string(ClusterMembership.ClusterID)==clusterID & ClusterMembership.Photoperiod_h==photo, :);
            if height(M) == 0, continue; end

            Y = [];
            yLabels = strings(0,1);
            uniqueFiles = unique(string(M.File), 'stable');

            for f = 1:numel(uniqueFiles)
                fileStem = uniqueFiles(f);
                rawPath = find_raw_path(RawFileIndex, fileStem);
                if strlength(rawPath) == 0
                    warnRows(end+1,:) = {clusterID, photo, fileStem, "", "Raw file not found"}; %#ok<AGROW>
                    log_line(LOG, 'WARNING: raw file not found for %s.', fileStem);
                    continue;
                end

                try
                    Traw = readtable(rawPath, 'VariableNamingRule','preserve');
                catch ME
                    warnRows(end+1,:) = {clusterID, photo, fileStem, "", "Could not read raw Excel file: " + string(ME.message)}; %#ok<AGROW>
                    continue;
                end

                [timeH, timeWarn] = extract_time_hours(Traw);
                if strlength(timeWarn) > 0
                    warnRows(end+1,:) = {clusterID, photo, fileStem, "", timeWarn}; %#ok<AGROW>
                end
                if isempty(timeH) || all(~isfinite(timeH))
                    warnRows(end+1,:) = {clusterID, photo, fileStem, "", "No usable time column"}; %#ok<AGROW>
                    continue;
                end

                Mf0 = M(string(M.File)==fileStem, :);
                if height(Mf0) == 0
                    continue;
                end

                % One activity trace per mouse/signal per cluster/photoperiod.
                % The activity-component filter is cluster-wide, so repeated
                % candidates from the same mouse would otherwise duplicate the
                % same trace.
                [~, uSignalIdx] = unique(string(Mf0.File) + "|" + string(Mf0.SignalID), 'stable');
                Mf = Mf0(uSignalIdx, :);

                for r = 1:height(Mf)
                    sig = string(Mf.SignalID(r));
                    col = find_table_column(Traw, sig);
                    if strlength(col) == 0
                        warnRows(end+1,:) = {clusterID, photo, fileStem, sig, "Signal/activity column not found in raw file"}; %#ok<AGROW>
                        continue;
                    end

                    rawY = to_double(Traw.(col));
                    [lowPad, highPad] = padded_period_window(double(C.PeriodLow_h), double(C.PeriodHigh_h), padFrac);
                    [comp, timeCompH, fWarn] = period_targeted_activity_component(timeH, rawY, lowPad, highPad, filterOrder, minCycles);
                    if strlength(fWarn) > 0
                        warnRows(end+1,:) = {clusterID, photo, fileStem, sig, fWarn}; %#ok<AGROW>
                        continue;
                    end

                    zt = mod(timeCompH(:), 24);
                    binned = NaN(1, numel(centers));
                    nobs = zeros(1, numel(centers));
                    for b = 1:numel(centers)
                        idx = zt >= edges(b) & zt < edges(b+1) & isfinite(comp);
                        if b == numel(centers)
                            idx = zt >= edges(b) & zt <= edges(b+1) & isfinite(comp);
                        end
                        binned(b) = mean(comp(idx), 'omitnan');
                        nobs(b) = sum(idx);
                        rows(end+1,:) = {clusterID, double(C.ClusterRank), string(C.BandName), photo, fileStem, sig, centers(b), edges(b), edges(b+1), ...
                            binned(b), nobs(b), lowPad, highPad}; %#ok<AGROW>
                    end
                    Y(end+1,:) = binned; %#ok<AGROW>
                    yLabels(end+1,1) = sig; %#ok<AGROW>
                end
            end

            if isempty(Y), continue; end

            f = figure('Color','w','Position',[100 100 950 520]);
            ax = axes(f); hold(ax,'on');
            for r = 1:size(Y,1)
                plot_background_line(ax, centers, Y(r,:), indivColour);
            end
            meanY = mean(Y, 1, 'omitnan');
            plot(ax, centers, meanY, '-', 'LineWidth', 2.8, 'Color', meanColour);
            yline(ax, 0, ':', 'Color', [0.3 0.3 0.3], 'HandleVisibility','off');
            add_light_dark_annotations(ax, photo, llValue, llRefPhoto, [0.82 0.82 0.82], [0.72 0.72 0.72], [0.05 0.05 0.05], 0.10);
            xlabel(ax, 'ZT / projected ZT (h)', 'FontWeight','bold');
            ylabel(ax, 'Activity (z-scored)', 'FontWeight','bold');
            title(ax, sprintf('24 h activity component | %s | %s | C%d %.2g-%.2g h', ...
                photo_label(photo), string(C.BandName), double(C.ClusterRank), double(C.PeriodLow_h), double(C.PeriodHigh_h)), ...
                'Interpreter','none', 'FontWeight','bold');
            xlim(ax, [0 24]);
            box(ax,'off');
            set(ax,'TickDir','out','FontName',figFont,'FontSize',14,'LineWidth',1.1);

            out = fullfile(figDir, sprintf('Activity_zscored_24h_%s_%s%s', photo_label(photo), sanitise_filename(clusterID), ext));
            exportgraphics(f, out, 'Resolution', dpi);
            close(f);
            PlotFiles = add_plotfile(PlotFiles, "Activity_zscored_24h", clusterID, photo_label(photo), out);
            log_line(LOG, 'Wrote %s', out);
        end
    end

    ActivityZTTable = rows_to_table(rows, hdr);
    Warnings = rows_to_table(warnRows, whdr);
end

function RawFileIndex = index_raw_excel_files(rawFolder)
    files = [dir(fullfile(rawFolder, '*.xlsx')); dir(fullfile(rawFolder, '*.xls'))];
    % Recursive fallback
    if isempty(files)
        files = [dir(fullfile(rawFolder, '**', '*.xlsx')); dir(fullfile(rawFolder, '**', '*.xls'))];
    else
        recFiles = [dir(fullfile(rawFolder, '**', '*.xlsx')); dir(fullfile(rawFolder, '**', '*.xls'))];
        files = [files; recFiles]; %#ok<AGROW>
    end

    rows = {};
    hdr = {'FileStem','FileName','Folder','FullPath'};
    seen = strings(0,1);
    for i = 1:numel(files)
        if files(i).isdir, continue; end
        [~, stem, ext] = fileparts(files(i).name);
        full = string(fullfile(files(i).folder, files(i).name));
        if any(seen == full), continue; end
        seen(end+1,1) = full; %#ok<AGROW>
        rows(end+1,:) = {string(stem), string(files(i).name), string(files(i).folder), full}; %#ok<AGROW>
    end
    RawFileIndex = rows_to_table(rows, hdr);
end

function rawPath = find_raw_path(RawFileIndex, fileStem)
    rawPath = "";
    if isempty(RawFileIndex) || height(RawFileIndex) == 0
        return;
    end

    fileStem = string(fileStem);
    stems = string(RawFileIndex.FileStem);

    idx = stems == fileStem;

    if ~any(idx)
        nFile = normalise_name(fileStem);
        nStem = strings(size(stems));
        for k = 1:numel(stems)
            nStem(k) = normalise_name(stems(k));
        end
        idx = nStem == nFile;
    end

    if ~any(idx)
        nFile = normalise_name(fileStem);
        nStem = strings(size(stems));
        for k = 1:numel(stems)
            nStem(k) = normalise_name(stems(k));
        end

        idx = false(size(stems));
        for k = 1:numel(stems)
            if strlength(nStem(k)) > 0
                idx(k) = contains(nStem(k), nFile) || contains(nFile, nStem(k));
            end
        end
    end

    if any(idx)
        cand = RawFileIndex(idx,:);
        rawPath = string(cand.FullPath(1));
    end
end

function [timeH, warn] = extract_time_hours(T)
    warn = "";
    timeH = [];

    names = string(T.Properties.VariableNames);
    lowerNames = lower(names);
    idx = find(contains(lowerNames, "time") | contains(lowerNames, "date"), 1, 'first');
    if isempty(idx)
        idx = 1;
        warn = "No explicit time column found; using first column as time.";
    end

    v = T.(names(idx));
    if isdatetime(v)
        timeH = hours(v - v(1));
    else
        x = to_double(v);
        nm = lowerNames(idx);
        if contains(nm, "day")
            timeH = x * 24;
        elseif contains(nm, "min")
            timeH = x / 60;
        else
            % Default for current mouse behaviour files: Time (hr)
            timeH = x;
        end
    end

    timeH = double(timeH(:));
    if any(diff(timeH(isfinite(timeH))) < 0)
        warn = warn + " Time column is not monotonic.";
    end
end

function col = find_table_column(T, signalID)
    col = "";
    names = string(T.Properties.VariableNames);
    signalID = string(signalID);

    exact = find(names == signalID, 1, 'first');
    if ~isempty(exact)
        col = names(exact);
        return;
    end

    normSignal = normalise_name(signalID);
    normNames = strings(size(names));
    for k = 1:numel(names)
        normNames(k) = normalise_name(names(k));
    end
    idx = find(normNames == normSignal, 1, 'first');
    if ~isempty(idx)
        col = names(idx);
        return;
    end

    validNames = strlength(normNames) > 0;
    idx = find(validNames & (contains(normNames, normSignal) | contains(normSignal, normNames)), 1, 'first');
    if ~isempty(idx)
        col = names(idx);
    end
end

function [compZ, timeOutH, warn] = period_targeted_activity_component(timeH, y, periodLowH, periodHighH, filterOrder, minCycles)
    warn = "";
    compZ = [];
    timeOutH = [];

    timeH = double(timeH(:));
    y = double(y(:));

    okTime = isfinite(timeH);
    if sum(okTime) < 20
        warn = "Too few finite time points for filtering.";
        return;
    end

    timeH = timeH(okTime);
    y = y(okTime);

    if sum(isfinite(y)) < 20
        warn = "Too few finite raw activity points for filtering.";
        return;
    end

    [timeH, ord] = sort(timeH);
    y = y(ord);
    timeOutH = timeH;

    dt = median(diff(timeH), 'omitnan');
    if ~isfinite(dt) || dt <= 0
        warn = "Could not estimate sampling interval.";
        return;
    end

    totalH = max(timeH) - min(timeH);
    if totalH < minCycles * periodHighH
        warn = sprintf('Recording too short for at least %d cycles of %.3g h period.', minCycles, periodHighH);
        return;
    end

    % Fill missing or irregular values on the observed time base only.
    y = fillmissing(y, 'pchip', 'MaxGap', max(3, round(2/dt)));
    y = fillmissing(y, 'linear');
    y = fillmissing(y, 'nearest');

    y = zscore_local(y);
    Fs = 1 / dt; % samples per hour
    nyq = Fs / 2;

    fLow = 1 / periodHighH;
    fHigh = 1 / periodLowH;

    if fHigh >= nyq
        warn = sprintf('Filter high frequency %.3g cycles/h exceeds Nyquist %.3g cycles/h.', fHigh, nyq);
        return;
    end
    if fLow <= 0 || fLow >= fHigh
        warn = "Invalid period-targeted filter frequency bounds.";
        return;
    end

    try
        [b,a] = butter(filterOrder, [fLow fHigh] / nyq, 'bandpass');
        if numel(y) <= 3 * max(numel(a), numel(b))
            warn = "Too few points for filtfilt with the requested bandpass filter.";
            return;
        end
        comp = filtfilt(b, a, y);
        compZ = zscore_local(comp);
    catch ME
        warn = "Bandpass filtering failed: " + string(ME.message);
        compZ = [];
    end
end

function add_light_dark_annotations(ax, photo, llValue, llRefPhoto, darkCol, projCol, lineCol, llAlpha)
    axes(ax); %#ok<LAXES>
    yl = ylim(ax);
    hold(ax,'on');

    isLL = photo >= llValue || photo >= 24;
    if isLL
        x1 = llRefPhoto;
        x2 = 24;
        ph = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], projCol, ...
            'FaceAlpha', llAlpha, 'EdgeColor', 'none', 'HandleVisibility','off');
        try, uistack(ph,'bottom'); catch, end
        xline(ax, 0, '-', 'Projected DL', 'Color', lineCol, 'LineWidth', 1.2, 'LabelOrientation','horizontal', 'HandleVisibility','off');
        xline(ax, llRefPhoto, '-', 'Projected LD', 'Color', lineCol, 'LineWidth', 1.2, 'LabelOrientation','horizontal', 'HandleVisibility','off');
    else
        if photo < 24
            ph = patch(ax, [photo 24 24 photo], [yl(1) yl(1) yl(2) yl(2)], darkCol, ...
                'FaceAlpha', 0.30, 'EdgeColor', 'none', 'HandleVisibility','off');
            try, uistack(ph,'bottom'); catch, end
        end
        xline(ax, 0, '-', 'DL', 'Color', lineCol, 'LineWidth', 1.2, 'LabelOrientation','horizontal', 'HandleVisibility','off');
        if photo > 0 && photo < 24
            xline(ax, photo, '-', 'LD', 'Color', lineCol, 'LineWidth', 1.2, 'LabelOrientation','horizontal', 'HandleVisibility','off');
        end
    end
    ylim(ax, yl);
end

function [preWin, postWin] = transition_window_limits(photoperiod_h, transitionType, periWindowH, adaptiveWindows, windowFraction)
    photo = double(photoperiod_h);
    dark_h = 24 - photo;
    preWin = periWindowH;
    postWin = periWindowH;

    if ~adaptiveWindows
        return;
    end

    transitionType = string(transitionType);
    switch transitionType
        case "DL"
            preCap  = windowFraction * dark_h;
            postCap = windowFraction * photo;
        case "LD"
            preCap  = windowFraction * photo;
            postCap = windowFraction * dark_h;
        otherwise
            preCap = periWindowH;
            postCap = periWindowH;
    end

    preWin  = min(periWindowH, max(0, preCap));
    postWin = min(periWindowH, max(0, postCap));
end

function [R, mu, csd] = circ_summary(theta)
    theta = theta(:);
    theta = theta(isfinite(theta));
    if isempty(theta)
        R = NaN; mu = NaN; csd = NaN;
        return;
    end
    z = mean(exp(1i*theta));
    R = abs(z);
    mu = angle(z);
    if R <= 0
        csd = Inf;
    else
        csd = sqrt(max(0, -2*log(R)));
    end
end

function plot_background_line(ax, x, y, rgb)
    try
        plot(ax, x, y, '-', 'Color', [rgb 0.35], 'LineWidth', 0.8, 'HandleVisibility','off');
    catch
        plot(ax, x, y, '-', 'Color', rgb, 'LineWidth', 0.8, 'HandleVisibility','off');
    end
end

function [low, high] = parse_band_bounds_h(band)
    band = string(band);
    nums = regexp(band, 'UR_(\d+)_?(\d+)?', 'tokens', 'once');
    if isempty(nums)
        low = NaN; high = NaN;
        return;
    end
    low = str2double(nums{1});
    high = str2double(nums{2});
end

function x = to_double(v)
    if isnumeric(v)
        x = double(v(:));
    elseif islogical(v)
        x = double(v(:));
    elseif iscell(v)
        x = str2double(string(v(:)));
    else
        x = str2double(string(v(:)));
    end
end

function x = to_logical(v)
    if islogical(v)
        x = v(:);
    elseif isnumeric(v)
        x = v(:) ~= 0;
    else
        s = lower(strtrim(string(v(:))));
        x = ismember(s, ["true","1","yes","y","pass","passed"]);
    end
end

function y = zscore_local(x)
    x = double(x(:));
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        y = x - mu;
    else
        y = (x - mu) ./ sd;
    end
end

function n = count_unique(x)
    if isempty(x)
        n = 0;
        return;
    end
    if isnumeric(x)
        x = x(isfinite(x));
        n = numel(unique(x));
    else
        x = string(x);
        x = x(~ismissing(x) & strlength(x) > 0);
        n = numel(unique(x));
    end
end

function s = robust_first_string(x)
    x = string(x);
    if isempty(x)
        s = "";
    else
        s = x(1);
    end
end

function out = normalise_name(s)
    out = lower(regexprep(string(s), '[^a-zA-Z0-9]', ''));
end

function q = local_iqr(x)
    x = x(:);
    x = x(isfinite(x));
    if isempty(x)
        q = NaN;
    else
        q = prctile(x,75) - prctile(x,25);
    end
end

function lab = photo_label(photo)
    if ~isfinite(photo)
        lab = 'LNA';
    elseif photo >= 24
        lab = 'LL';
    else
        lab = sprintf('L%.3g', photo);
    end
end

function s = safe_period(x)
    s = regexprep(sprintf('%.2f', x), '\.', 'p');
end

function s = sanitise_filename(s)
    s = string(s);
    s = regexprep(s, '[^\w\-.]+', '_');
    s = regexprep(s, '_+', '_');
    s = char(s);
end

function nm = transition_display_name(t)
    t = string(t);
    switch t
        case "ProjectedDL_LL"
            nm = "Projected DL";
        case "ProjectedLD_LL"
            nm = "Projected LD";
        otherwise
            nm = t;
    end
end

function T = rows_to_table(rows, hdr)
    % Convert a cell-row accumulator into a table with sensible variable
    % types. This avoids later failures where numeric columns remain as cell
    % arrays after cell2table.
    if isempty(rows)
        T = cell2table(cell(0, numel(hdr)), 'VariableNames', hdr);
        return;
    end

    T = cell2table(rows, 'VariableNames', hdr);

    for j = 1:numel(hdr)
        col = rows(:, j);

        isNumLike = true(size(col));
        for k = 1:numel(col)
            x = col{k};
            isNumLike(k) = isempty(x) || ...
                (isnumeric(x) && isscalar(x)) || ...
                (islogical(x) && isscalar(x));
        end

        if all(isNumLike)
            vals = NaN(numel(col), 1);
            for k = 1:numel(col)
                x = col{k};
                if isempty(x)
                    vals(k) = NaN;
                else
                    vals(k) = double(x);
                end
            end
            T.(hdr{j}) = vals;
            continue;
        end

        isTextLike = true(size(col));
        for k = 1:numel(col)
            x = col{k};
            isTextLike(k) = isempty(x) || ischar(x) || isstring(x) || iscategorical(x);
        end

        if all(isTextLike)
            vals = strings(numel(col), 1);
            for k = 1:numel(col)
                x = col{k};
                if isempty(x)
                    vals(k) = "";
                else
                    vals(k) = string(x);
                end
            end
            T.(hdr{j}) = vals;
        end
    end
end

function PF = empty_plotfiles()
    PF = table(string.empty(0,1), string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
        'VariableNames', {'PlotType','ClusterID','Photoperiod','FilePath'});
end

function PF = add_plotfile(PF, plotType, clusterID, photoperiod, filePath)
    newRow = cell2table({string(plotType), string(clusterID), string(photoperiod), string(filePath)}, ...
        'VariableNames', PF.Properties.VariableNames);
    PF = [PF; newRow];
end

function safe_writetable(T, xlsx, sheet)
    if isempty(T)
        T = table();
    end
    sheet = char(string(sheet));
    if strlength(string(sheet)) > 31
        sheet = char(extractBefore(string(sheet), 32));
    end
    try
        writetable(T, xlsx, 'Sheet', sheet);
    catch ME
        warning('Could not write sheet %s: %s', sheet, ME.message);
    end
end

function safe_delete_file(p)
    if isfile(p)
        try
            delete(p);
        catch
        end
    end
end

function ensure_dir(d)
    if ~exist(d, 'dir')
        mkdir(d);
    end
end

function log_line(fid, varargin)
    if fid > 0
        fprintf(fid, [varargin{1} '\n'], varargin{2:end});
    end
end

function fclose_if_valid(fid)
    if ~isempty(fid) && fid > 0
        try
            fclose(fid);
        catch
        end
    end
end
