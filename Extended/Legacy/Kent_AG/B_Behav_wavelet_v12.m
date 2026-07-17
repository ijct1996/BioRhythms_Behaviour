%% ========================================================================
% Script 1 (v12 - validated Raw ultradian pipeline + ridge phase handoff)
% File: Behav_wavelet_v12.m
% ========================================================================
% Purpose
%   Per input file:
%     - Create per-file output folder under parent output
%     - Wavelet analysis (Raw / HSub / Raw+HSub)
%     - Light/Dark + Phase derived ONLY from Light duration (h), ZT0=lights-on
%     - Outputs: per-file Excel + figures + logs
%     - Handoff outputs (single shared handoff folder):
%         1) WP_Summary__<fileStem>.mat  (small summary tables + meta)
%         2) WP_TS__<fileStem>.mat       (time-series long tables only, with ValidFlag)
%
% v5.6 simplified additions (only):
%   1) Conditions-by-columns UI in preflight, plus ConditionMap workbook
%   2) Per-file averaged Abs scalograms across signals (per condition + overall), Raw/HSub
%   3) Across-photoperiod pooled averaged Abs scalograms when multiple files selected,
%      grouped by photoperiod and overall, Raw/HSub, per condition + overall
%
% v5.7 REQUIRED CHANGES (STRICT):
%   A) Across-photoperiod averaged Abs scalograms:
%      - Replace per photoperiod separate plots with a single combined scalogram per group.
%      - Combined scalogram concatenates photoperiod blocks along the time axis.
%      - Mark each photoperiod boundary with a white dotted vertical line.
%      - Apply to Raw and HSub.
%      - Save only JPEGs (no MAT files).
%   B) Raw vs HSub intensity scaling (averaged Abs scalograms only):
%      - For per-file and across-photoperiod averaged Abs scalograms:
%        enforce identical CLim scaling between Raw and HSub for the same group.
%      - CLimMax = robust 99th percentile across BOTH Raw and HSub averages for that group.
%      - Apply caxis([0 CLimMax]) to both.
%   C) Remove averaged-scalogram MAT outputs:
%      - Do not save any MAT files for per-file or cross-file averaged Abs scalograms.
%      - Only save JPEG figures.
%
% v5.8 ADDITION:
%   - Use Paul Tol (Bright) colour palette consistently for band plots and ridge overlays.
%   - Write a sidecar band colour map workbook into the HANDOFF folder (no pkgS/pkgTS changes):
%       AcrossPhotoperiod_Input/WP_BandColourMap.xlsx
%
% v12 ADDITIONS FOR VALIDATED-RAW ULTRADIAN PIPELINE:
%   - Raw data are analysed for BOTH CR and UR bands.
%   - HSub is used as a validation layer, not as the final biological signal.
%   - HSub selection now prefers the v12 ValidationManifest sheet and uses
%     PrimaryValidationResidual = SEL_P360 when available.
%   - Adds PeriodCandidates_Long to the summary handoff.
%   - Adds RidgeTrajectory_Long and RidgePhase_Long to the TS handoff.
%   - Extracts ridge-following wavelet phase from the complex CWT coefficient.
%
% HARD CONSTRAINTS:
%   - Retain all existing functionalities and outputs EXCEPT where explicitly requested above.
%   - Do NOT remove or rename existing pkgS/pkgTS structures, field names, or long-table schemas consumed by Script 2.
%   - New tables are additive for validation and resynchronisation scripts.
%   - Do NOT change wavelet period limits, bands, ridge extraction logic, or phase split logic except Raw now also exports UR bands.
%   - All output figures must remain JPEG at SAVE_DPI = 600.
%
% MATLAB: R2025a/R2025b
% ========================================================================

clearvars; close all; clc;

%% ----------------------------- NAMING -----------------------------------
NAMES = struct();

% Handoff folder (for Script 2)
NAMES.HANDOFF_DIR       = 'AcrossPhotoperiod_Input';
NAMES.SUM_PREFIX        = 'WP_Summary__%s.mat';
NAMES.TS_PREFIX         = 'WP_TS__%s.mat';
NAMES.HANDOFF_INDEX     = 'WP_Package_Index.xlsx';

% Per-file Excel outputs
NAMES.XLSX_SUMMARY      = 'WP_Summary.xlsx';
NAMES.XLSX_DETAIL       = 'WP_Detail.xlsx';
NAMES.XLSX_BANDTS       = 'WP_BandTimeSeries.xlsx';

% Per-file folder names
NAMES.REPORTS_DIR       = 'Reports';
NAMES.FIG_DIR           = 'Figures';
NAMES.LOGS_DIR          = 'Logs';
NAMES.WIDE_DIR          = 'WideTimeSeries';
NAMES.DETAILCACHE_DIR   = 'DetailCache';
NAMES.PLOTMETA_XLSX     = 'WP_PlotMeta.xlsx';

% Conditions workbook
NAMES.CONDITIONS_DIR    = 'Conditions';
NAMES.CONDITIONMAP_XLSX = 'WP_ConditionMap.xlsx';

% Wavelet scalograms subfolders
NAMES.FIG_WAVELET       = 'Wavelet';
NAMES.FIG_WAV_ABS       = 'Abs';
NAMES.FIG_WAV_NORM      = 'NormZ';
NAMES.FIG_WAV_BANDS     = 'BandLines';
NAMES.FIG_WAV_RIDGE     = 'RidgeOverlay';

% Averaged scalograms (per file)
NAMES.FIG_WAV_AVG       = 'Averages';

% Spectra subfolders
NAMES.FIG_POWER         = 'Power';
NAMES.PWR_RAW           = 'Raw';
NAMES.PWR_HSUB          = 'HSub';
NAMES.PWR_GLOBAL        = 'Global';
NAMES.PWR_COND          = 'Cond';

% Band/Ridge time-series subfolders
NAMES.FIG_BANDRIDGE     = 'BandRidge';
NAMES.BR_BP             = 'BandPower';
NAMES.BR_RP             = 'RidgePeriod';
NAMES.BR_RPOW           = 'RidgePower';
NAMES.BR_COMBINED       = 'Combined';

% Figure filenames
NAMES.FN_WAV_ABS        = 'WAV_%s_%s_%s_Abs.jpg';
NAMES.FN_WAV_NORM       = 'WAV_%s_%s_%s_NormZ.jpg';
NAMES.FN_WAV_BANDS      = 'WAV_%s_%s_%s_Bands.jpg';
NAMES.FN_WAV_RIDGE      = 'WAV_%s_%s_%s_Ridge.jpg';

NAMES.FN_PWR_GLOB       = 'PWR_%s_%s_%s.jpg';
NAMES.FN_PWR_COND       = 'PWR_%s_%s_%s_Cond_%s.jpg';

NAMES.FN_BP_COMB        = 'BP_%s_%s_%s.jpg';
NAMES.FN_RP_COMB        = 'RP_%s_%s_%s.jpg';
NAMES.FN_RPOW_COMB      = 'RPOW_%s_%s_%s.jpg';

% Per-file averaged scalogram previews
NAMES.FN_AVG_ABS        = 'AVGABS_%s_%s_%s.jpg'; % (FileStem, Source, GroupLabel)

SAVE_DPI = 600;

%% ----------------------------- CONSTANTS --------------------------------
PERIOD_MIN_MIN = 60;
PERIOD_MAX_MIN = 1590;                 % 26.5 h max (DO NOT CHANGE)
GLOBAL_YMAX_HR = hours(minutes(PERIOD_MAX_MIN));
SCALO_YTICKS_H = 0:4:26;
SCALO_HLINES_H = [1 3 6 9 12 18 20 28];

TOPN = 5;

BAND_CR = [20, 28];
BANDS_UR = [
    1, 3;
    3, 6;
    6, 9;
    9,12;
    12,18
];

% v12: Raw is now the biological source for CR and UR. HSub is validation only.
BANDS_RAW_ALL = [reshape(BAND_CR,1,2); BANDS_UR];
BAND_NAMES_ALL = {'CR_20_28','UR_1_3','UR_3_6','UR_6_9','UR_9_12','UR_12_18'};
BAND_NAMES_UR  = {'UR_1_3','UR_3_6','UR_6_9','UR_9_12','UR_12_18'};

% v12 candidate/QC defaults. These are exported in README/meta and can be changed centrally.
HS_PRIMARY_VALIDATION_RESIDUAL   = 'SEL_P360';
HS_SECONDARY_VALIDATION_RESIDUAL = 'SEL_P60';
CANDIDATE_MIN_RIDGE_COVERAGE    = 0.50;
CANDIDATE_MIN_COI_VALID_FRAC     = 0.50;

FORCE_ZERO_Y_BP_RPOW = true;
FORCE_ZERO_Y_RP      = true;

% Cross-file pooled averages (plain Abs only)
DO_CROSSFILE_POOLED_AVG = true;

%% ----------------------------- COLOURS ----------------------------------
% Paul Tol (Bright) qualitative scheme (colour-blind safe)
% Source: https://sronpersonalpages.nl/~pault/ (Bright)
TOL = struct();
TOL.BrightHex = {'#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB'};
TOL.BrightRGB = hex2rgb_array(TOL.BrightHex);

% Band colour map (single source of truth, stable names)
BAND_COLOUR = struct();
BAND_COLOUR.CR_20_28 = TOL.BrightRGB(1,:);  % blue
BAND_COLOUR.UR_1_3   = TOL.BrightRGB(2,:);  % red/pink
BAND_COLOUR.UR_3_6   = TOL.BrightRGB(3,:);  % green
BAND_COLOUR.UR_6_9   = TOL.BrightRGB(4,:);  % yellow
BAND_COLOUR.UR_9_12  = TOL.BrightRGB(5,:);  % cyan
BAND_COLOUR.UR_12_18 = TOL.BrightRGB(6,:);  % purple

%% ========================================================================
% 1) SELECT INPUT FILE(S)
% ========================================================================
runMode = questdlg( ...
    'Choose input mode:', ...
    'Wavelet analysis: input mode', ...
    'Single file', 'Multiple files', 'Cancel', ...
    'Single file');

if isempty(runMode) || strcmpi(runMode,'Cancel')
    fprintf('No mode selected. Exiting.\n');
    return;
end

[fileList, ~] = select_input_files(runMode);
if isempty(fileList)
    fprintf('No input file(s) selected. Exiting.\n');
    return;
end

nFiles = numel(fileList);
fprintf('Selected %d input file(s).\n', nFiles);

%% ========================================================================
% 2) RUN TYPE
% ========================================================================
runType = questdlg( ...
    'Is this run Raw only, HSub only, or Raw + HSub?', ...
    'Wavelet analysis: run type', ...
    'Raw only', 'HSub only', 'Raw + HSub', ...
    'Raw + HSub');

if isempty(runType)
    fprintf('No run type selected. Exiting.\n');
    return;
end

doRaw  = strcmpi(runType, 'Raw only') || strcmpi(runType, 'Raw + HSub');
doHSub = strcmpi(runType, 'HSub only') || strcmpi(runType, 'Raw + HSub');

%% ========================================================================
% 3) HSub MODE (Auto vs Manual)
% ========================================================================
HS_MODE = 'Auto';
HS_PARENT = '';
MANUAL = struct();
MANUAL.summaryFile   = '';
MANUAL.residualRoot  = '';
MANUAL.residualFiles = {};
MANUAL.map           = struct();
MANUAL.mapOK         = false;

HSIndex = struct();
HSIndex.ok = false;
HSIndex.byStem = struct();

if doHSub
    hsMode = questdlg( ...
        'HSub selection mode:', ...
        'HSub mode', ...
        'Auto (select HSub parent folder)', ...
        'Manual (select HS summary + residual workbooks)', ...
        'Auto (select HSub parent folder)', ...
        'Auto (select HSub parent folder)');

    if isempty(hsMode)
        fprintf('No HSub mode selected. Falling back to Raw only.\n');
        doHSub = false;
    else
        if startsWith(lower(hsMode),'auto')
            HS_MODE = 'Auto';
        else
            HS_MODE = 'Manual';
        end
    end
end

%% ========================================================================
% 4) PHASE SPLIT (Light vs Dark)  ZT0 = lights-on (fact)
% ========================================================================
ansPhase = questdlg( ...
    'Enable phase split (Light vs Dark) derived from Light duration (h) and ZT0=lights-on?', ...
    'Phase split', ...
    'Yes', 'No', 'No');

DO_PHASE_SPLIT = strcmpi(ansPhase, 'Yes');

%% ========================================================================
% 5) SELECT PARENT OUTPUT FOLDER
% ========================================================================
parentOut = uigetdir(pwd, 'Select PARENT output folder (per-file subfolder created for each input file)');
if isequal(parentOut,0)
    fprintf('No output folder selected. Exiting.\n');
    return;
end

% Shared handoff folder
handoffDir = fullfile(parentOut, NAMES.HANDOFF_DIR);
ensure_dir(handoffDir);

% Write sidecar band colour map (does NOT touch pkgS/pkgTS)
bandColourMapXLSX = fullfile(handoffDir, 'WP_BandColourMap.xlsx');
write_band_colour_map_xlsx(bandColourMapXLSX, BAND_COLOUR, TOL);

handoffIndexRows = {};
handoffIndexHdr  = {'FileStem','InputFile','PerFileFolder','SummaryMat','TSMat','RunType','HSubMode','PhaseSplitEnabled'};

% Cross-file averages root (additive)
crossAvgRoot = fullfile(parentOut, 'AcrossFile_Averages');
ensure_dir(crossAvgRoot);

% Combined across-photoperiod output root (required)
crossCombinedRoot = fullfile(crossAvgRoot, 'Combined');
ensure_dir(crossCombinedRoot);

%% ========================================================================
% 6) HSub setup (once)
% ========================================================================
if doHSub
    if strcmpi(HS_MODE,'Auto')
        HS_PARENT = uigetdir(pwd, 'Select HSub parent folder (e.g. "01. HSub")');
        if isequal(HS_PARENT,0)
            fprintf('No HSub parent folder selected. Falling back to Raw only.\n');
            doHSub = false;
        else
            fprintf('\nAuto HSub: indexing HS_Summary.xlsx under:\n  %s\n', HS_PARENT);
            HSIndex = build_hsub_index_from_parent(HS_PARENT);
            if ~HSIndex.ok
                fprintf('Auto HSub indexing failed. Falling back to Raw only.\n');
                doHSub = false;
            end
        end
    else
        fprintf('\nManual HSub: select residual time-series workbook(s) (typically 4)...\n');
        [resFiles, resPath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, ...
            'Manual HSub: Select residual time-series workbook(s) (typically 4)', ...
            'MultiSelect','on');

        if isequal(resFiles,0)
            fprintf('No residual workbooks selected. Falling back to Raw only.\n');
            doHSub = false;
        else
            if ischar(resFiles), resFiles = {resFiles}; end
            resFiles = resFiles(:);
            MANUAL.residualFiles = cell(numel(resFiles),1);
            for i = 1:numel(resFiles)
                MANUAL.residualFiles{i} = fullfile(resPath, resFiles{i});
            end

            hsRoot = uigetdir(resPath, 'Manual optional: select residual ROOT folder (parent containing TS\RES...). Cancel to skip.');
            if ~isequal(hsRoot,0)
                MANUAL.residualRoot = hsRoot;
            else
                MANUAL.residualRoot = '';
            end

            fprintf('\nManual HSub: select HS summary workbook...\n');
            [hsSummaryName, hsSummaryPath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, ...
                'Manual HSub: Select HS summary workbook (HS_Summary.xlsx)');

            if isequal(hsSummaryName,0) || isequal(hsSummaryPath,0)
                fprintf('No HS summary selected. Falling back to Raw only.\n');
                doHSub = false;
            else
                MANUAL.summaryFile = fullfile(hsSummaryPath, hsSummaryName);
                [MANUAL.map, MANUAL.mapOK] = load_hsub_recommend_map(MANUAL.summaryFile);
                if ~MANUAL.mapOK
                    fprintf('Recommend map could not be loaded. Falling back to Raw only.\n');
                    doHSub = false;
                end
            end
        end
    end
end

%% ========================================================================
% 7) PREFLIGHT mapping (Time + Light duration + Data + Conditions-by-columns)
% ========================================================================
jobs = repmat(struct(), nFiles, 1);

for iF = 1:nFiles
    inputFile = fileList{iF};
    [~, fileStem, ext] = fileparts(inputFile);
    fileBaseWithExt = [fileStem ext];

    fprintf('\nPreflight %d/%d: %s\n', iF, nFiles, fileBaseWithExt);

    [Tpf, ~] = read_input_table_preserve_robust(inputFile);
    [Tpf, ~] = drop_empty_columns_robust(Tpf);

    if isempty(Tpf) || width(Tpf) < 2
        error('Insufficient columns after import cleanup in "%s".', fileBaseWithExt);
    end

    mapPF = preflight_column_mapping_dialog_time_lightduration_and_conditions(Tpf, iF, nFiles);

    jobs(iF).inputFile = inputFile;
    jobs(iF).fileStem  = fileStem;
    jobs(iF).fileBaseWithExt = fileBaseWithExt;
    jobs(iF).preflightMap = mapPF;

    jobs(iF).outDir = fullfile(parentOut, sprintf('File_%02d_%s', iF, sanitise_filename(fileStem)));
end

%% ========================================================================
% 8) PROCESS EACH FILE
% ========================================================================
residCache = containers.Map('KeyType','char','ValueType','any');

% Cross-file pooled accumulators (Abs only), for Raw and HSub
CrossPool = init_crosspool();

for iF = 1:nFiles
    job = jobs(iF);

    fprintf('\n============================================================\n');
    fprintf('Processing file %d/%d: %s\n', iF, nFiles, job.fileBaseWithExt);
    fprintf('Per-file output folder:\n  %s\n', job.outDir);
    fprintf('============================================================\n');

    ensure_dir(job.outDir);

    dirs = build_perfile_dirs(job.outDir, NAMES);
    ensure_perfile_dirs(dirs);

    runLogPath = fullfile(dirs.Logs, 'LOG_Run.txt');
    if exist(runLogPath,'file'), delete(runLogPath); end
    RUNLOG = log_open(runLogPath);

    log_line(RUNLOG, 'Wavelet_analysis_v5_8(simplified) started: %s', datestr(now,31));
    log_line(RUNLOG, 'InputFile: %s', job.inputFile);
    log_line(RUNLOG, 'RunType: %s | doRaw=%d | doHSub=%d', runType, doRaw, doHSub);
    log_line(RUNLOG, 'PhaseSplitEnabled: %d', DO_PHASE_SPLIT);
    log_line(RUNLOG, 'HSubMode: %s', HS_MODE);

    % Resolve HS context for this file
    HS = struct();
    HS.enabled = doHSub;
    HS.summaryFile = '';
    HS.residualRoot = '';
    HS.residualFiles = {};
    HS.map = struct();
    HS.mapOK = false;

    if doHSub
        if strcmpi(HS_MODE,'Auto')
            [entry, found] = hsub_index_lookup(HSIndex, job.fileStem);
            if ~found
                HS.enabled = false;
                log_line(RUNLOG, 'HSub AUTO: no mapping found for stem %s. HSub disabled for this file.', job.fileStem);
            else
                HS.summaryFile   = entry.hsSummaryPath;
                HS.residualRoot  = entry.runFolder;
                HS.residualFiles = entry.residualFiles;
                [HS.map, HS.mapOK] = load_hsub_recommend_map(HS.summaryFile);
                if ~HS.mapOK
                    HS.enabled = false;
                    log_line(RUNLOG, 'HSub AUTO: Recommend map not OK. HSub disabled for this file.');
                end
            end
        else
            if MANUAL.mapOK
                HS.summaryFile   = MANUAL.summaryFile;
                HS.residualRoot  = MANUAL.residualRoot;
                HS.residualFiles = MANUAL.residualFiles;
                HS.map           = MANUAL.map;
                HS.mapOK         = true;
            else
                HS.enabled = false;
            end
        end
    end

    % Per-file Excel outputs
    summaryXLSX   = fullfile(dirs.Reports, NAMES.XLSX_SUMMARY);
    detailXLSX    = fullfile(dirs.Reports, NAMES.XLSX_DETAIL);
    bandTSXLSX    = fullfile(dirs.Reports, NAMES.XLSX_BANDTS);
    plotMetaXLSX  = fullfile(dirs.Reports, NAMES.PLOTMETA_XLSX);

    delete_if_exists(summaryXLSX);
    delete_if_exists(detailXLSX);
    delete_if_exists(bandTSXLSX);
    delete_if_exists(plotMetaXLSX);

    detailCacheDir = dirs.DetailCache;
    ensure_dir(detailCacheDir);

    % ConditionMap workbook
    condDir = fullfile(dirs.Reports, NAMES.CONDITIONS_DIR);
    ensure_dir(condDir);
    condMapXLSX = fullfile(condDir, NAMES.CONDITIONMAP_XLSX);
    delete_if_exists(condMapXLSX);

    % Per-file accumulators (existing)
    peakRows = {};
    peakHdr  = {'File','SignalID','Source','ConditionParsed','Photoperiod_h','LightStateValue','Phase', ...
                'PeakRank','PeakPeriod_hr','PeakValue_log10','PeakProminence','PeakWidth'};

    inputsRows = {};
    inputsHdr  = {'File','TimeColumn','LightDurationColumn','Photoperiod_h', ...
                  'NumDataCols','ExcludedCols','Notes','RunType','HS_Summary','ResidualRoot','ResidualFiles', ...
                  'PhaseSplitEnabled','HSubMode','HSubEnabledForFile'};

    bandRows = {};
    bandHdr = {'File','SignalID','Source','Photoperiod_h', ...
               'MeanLogPower_CR', ...
               'MeanLogPower_UR_1_3','MeanLogPower_UR_3_6','MeanLogPower_UR_6_9','MeanLogPower_UR_9_12','MeanLogPower_UR_12_18'};

    bandPowRows  = {};
    ridgePRows   = {};
    ridgePowRows = {};
    bandTSHdr = {'File','SignalID','Source','BandName','Time_days','ZT_hours','LightStateValue','Phase','Value','ValidFlag'};

    bandCondRows = {};
    bandCondHdr = {'File','SignalID','Source','Photoperiod_h','LightStateValue','Phase','BandName', ...
                   'MeanBandPower_log10','SDBandPower_log10', ...
                   'MeanBandPower_linear','SDBandPower_linear', ...
                   'FracBandPower_linear', ...
                   'MeanRidgePeriod','SDRidgePeriod', ...
                   'MeanRidgePower_log10','SDRidgePower_log10'};

    % v12 validation/resynchronisation handoff tables. These are additive and do not
    % replace the legacy long-table outputs consumed by Script 2.
    periodCandidateRows = {};
    periodCandidateHdr = {'File','SignalID','ConditionParsed','Source','HSubResidualMode','Photoperiod_h','Phase','BandName', ...
                          'CandidateID','CandidateRank','MedianRidgePeriod_h','IQR_RidgePeriod_h', ...
                          'MeanBandPower_log10','SDBandPower_log10','MeanRidgePower_log10','SDRidgePower_log10', ...
                          'RidgeCoverageFrac','COIValidFrac','ValidPointCount','TotalPointCount','PassQC','QCReason'};

    ridgeTrajectoryRows = {};
    ridgeTrajectoryHdr = {'File','SignalID','ConditionParsed','Source','HSubResidualMode','Photoperiod_h','BandName','CandidateID', ...
                          'Time_days','ZT_hours','LightStateValue','Phase','BandPower_log10','RidgePeriod_h','RidgePower_log10','ValidFlag'};

    ridgePhaseRows = {};
    ridgePhaseHdr = {'File','SignalID','ConditionParsed','Source','HSubResidualMode','Photoperiod_h','BandName','CandidateID', ...
                     'Time_days','ZT_hours','LightStateValue','Phase','RidgePeriod_h','RidgePower_log10','RidgePhase_rad','ValidFlag'};

    plotMetaRows = {};
    plotMetaHdr = {'File','SignalID','Source','ScalogramType','CLimMin','CLimMax','PowerMax_abs','PowerP99_abs','Notes'};

    % Read input
    [T, ~] = read_input_table_preserve_robust(job.inputFile);
    [T, droppedNow] = drop_empty_columns_robust(T);
    if isempty(T) || width(T) < 2
        log_line(RUNLOG, 'SKIP: empty/insufficient columns after cleanup.');
        log_close(RUNLOG);
        continue;
    end

    [timeIdx, lightDurIdx, dataIdx, excludedNames, notes, condMap] = resolve_preflight_mapping_time_lightdur_and_conditions(T, job.preflightMap);

    varNames = T.Properties.VariableNames;
    timeName = varNames{timeIdx};
    lightDurName = varNames{lightDurIdx};

    if isempty(notes), notes = {}; end
    if ~isempty(droppedNow.DroppedNames)
        notes{end+1} = ['Dropped empty columns: ' strjoin(droppedNow.DroppedNames, ', ')];
    end

    % Write ConditionMap workbook
    write_condition_map_workbook(condMapXLSX, condMap, job.fileStem);

    timeCol = T{:, timeIdx};
    [timeMinutesAll, TsMinutes] = infer_time_minutes(timeCol, timeName);
    if ~isfinite(TsMinutes) || TsMinutes <= 0
        log_line(RUNLOG, 'SKIP: invalid sampling interval.');
        log_close(RUNLOG);
        continue;
    end

    N = height(T);
    time_min = timeMinutesAll(:);
    if numel(time_min) ~= N
        N = min(N, numel(time_min));
        time_min = time_min(1:N);
    end

    time_day = time_min / (60*24);
    time_hr  = time_min / 60;
    ZT_hr    = mod(time_hr,24);

    lightDurVec = to_numeric_vector(T{1:N, lightDurIdx});
    lightDurVec = standardise_signal_to_length(lightDurVec, N);
    photoperiod_h = robust_scalar_mode(lightDurVec);

    [isLight, lightStateStrVec, lightStateNumVec] = derive_light_state_from_lightdur(time_hr, lightDurVec);
    phaseMasks = struct();
    phaseMasks.All = true(N,1);
    if DO_PHASE_SPLIT
        phaseMasks.Light = isLight(:);
        phaseMasks.Dark  = ~phaseMasks.Light;
    end

    condChangeIdx = compute_changes_from_lightduration(lightDurVec);

    inputsRows(end+1,:) = { ...
        job.fileStem, timeName, lightDurName, photoperiod_h, ...
        numel(dataIdx), strjoin(excludedNames, ', '), strjoin(notes, ' | '), ...
        runType, HS.summaryFile, HS.residualRoot, strjoin(HS.residualFiles,' | '), ...
        DO_PHASE_SPLIT, HS_MODE, double(HS.enabled) ...
        }; %#ok<SAGROW>

    if ~(exist('cwtfilterbank','file')==2 && exist('cwt','file')==2)
        error('Wavelet Toolbox functions cwtfilterbank/cwt not found. Cannot proceed.');
    end

    FB = cwtfilterbank('SignalLength', N, ...
        'SamplingPeriod', minutes(TsMinutes), ...
        'PeriodLimits', [minutes(PERIOD_MIN_MIN), minutes(PERIOD_MAX_MIN)], ...
        'Wavelet', 'amor');

    SPEC_TEMPLATE = make_spec_template();

    % Wide TS workbooks (per file)
    wideBP_Raw   = table(time_day(:), ZT_hr(:), string(lightStateStrVec(:)), 'VariableNames', {'Time_days','ZT_hours','LightStateValue'});
    wideRP_Raw   = table(time_day(:), ZT_hr(:), string(lightStateStrVec(:)), 'VariableNames', {'Time_days','ZT_hours','LightStateValue'});
    wideRPOW_Raw = table(time_day(:), ZT_hr(:), string(lightStateStrVec(:)), 'VariableNames', {'Time_days','ZT_hours','LightStateValue'});

    wideBP_HSub   = table(time_day(:), ZT_hr(:), string(lightStateStrVec(:)), 'VariableNames', {'Time_days','ZT_hours','LightStateValue'});
    wideRP_HSub   = table(time_day(:), ZT_hr(:), string(lightStateStrVec(:)), 'VariableNames', {'Time_days','ZT_hours','LightStateValue'});
    wideRPOW_HSub = table(time_day(:), ZT_hr(:), string(lightStateStrVec(:)), 'VariableNames', {'Time_days','ZT_hours','LightStateValue'});

    % Detail spectra cache index
    detailCounter = 0;
    detailIndexRows = {};
    detailIndexHdr = {'ShortID','File','SignalID','Source','SpectrumSheet'};

    % For spectra plotting: keep per-column template store
    streamStoreRaw  = repmat(SPEC_TEMPLATE, numel(dataIdx), 1);
    streamStoreHSub = repmat(SPEC_TEMPLATE, numel(dataIdx), 1);
    globalPowerMax_Raw  = -inf;
    globalPowerMax_HSub = -inf;

    %% ----------------------- FIRST PASS PER COLUMN -----------------------
    for c = 1:numel(dataIdx)
        colIdx = dataIdx(c);
        signalID = char(string(varNames{colIdx}));

        log_line(RUNLOG, 'Column %d/%d: %s', c, numel(dataIdx), signalID);

        rec = struct(); rec.found = false;
        if HS.enabled && HS.mapOK
            rec = get_recommendation_for_column(HS.map, signalID);
            if ~rec.found
                log_line(RUNLOG, '  WARNING: No Recommend entry for "%s". HSub will be skipped for this mouse.', signalID);
            end
        end

        % ---------------- RAW ----------------
        if doRaw
            xRaw = standardise_signal_to_length(to_numeric_vector(T{1:N, colIdx}), N);
            if ~isempty(xRaw)
                xRaw(~isfinite(xRaw)) = 0;

                [specOut, peaksOut, plotMetaAdd, detailCounter, detailIndexRows] = process_one_stream( ...
                    xRaw, FB, N, time_day, condChangeIdx, ...
                    job.fileStem, signalID, 'Raw', dirs, NAMES, SAVE_DPI, ...
                    lightStateNumVec, photoperiod_h, ...
                    SPEC_TEMPLATE, SCALO_YTICKS_H, SCALO_HLINES_H, ...
                    BAND_CR, BANDS_UR, detailCacheDir, detailCounter, detailIndexRows, ...
                    BAND_COLOUR);

                if ~isempty(plotMetaAdd), plotMetaRows = [plotMetaRows; plotMetaAdd]; end %#ok<AGROW>

                streamStoreRaw(c) = specOut;
                if specOut.ok
                    globalPowerMax_Raw = max(globalPowerMax_Raw, specOut.localPowerMax);
                end
                if ~isempty(peaksOut), peakRows = [peakRows; peaksOut]; end %#ok<AGROW>

                if specOut.ok
                    % band summaries
                    bandRow = compute_band_summary_row(job.fileStem, signalID, 'Raw', photoperiod_h, specOut.periods_hours, specOut.avgPowerSpectrum, BAND_CR, BANDS_UR);
                    if ~isempty(bandRow), bandRows(end+1,:) = bandRow; end %#ok<SAGROW>

                    % v12: Raw exports CR + UR bands. Raw remains the biological signal.
                    bandNames = BAND_NAMES_ALL;
                    if isfield(specOut,'defaultBandTS') && ~isempty(specOut.defaultBandTS)
                        bandTS = specOut.defaultBandTS;
                    else
                        bandTS = compute_band_timeseries(specOut.periods_hours, specOut.logPow, BANDS_RAW_ALL, specOut.coi_hours);
                    end

                    [bandPowRows, ridgePRows, ridgePowRows] = append_band_timeseries_rows_with_valid( ...
                        bandPowRows, ridgePRows, ridgePowRows, ...
                        job.fileStem, signalID, 'Raw', bandNames, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTS);

                    bandCondRows = [bandCondRows; append_band_condition_summary_linearfrac( ...
                        job.fileStem, signalID, 'Raw', photoperiod_h, phaseMasks, bandNames, bandTS)]; %#ok<AGROW>

                    periodCandidateRows = [periodCandidateRows; append_period_candidate_rows( ...
                        job.fileStem, signalID, 'Raw', 'NA', photoperiod_h, phaseMasks, bandNames, bandTS, ...
                        CANDIDATE_MIN_RIDGE_COVERAGE, CANDIDATE_MIN_COI_VALID_FRAC)]; %#ok<AGROW>

                    ridgeTrajectoryRows = [ridgeTrajectoryRows; append_ridge_trajectory_rows( ...
                        job.fileStem, signalID, 'Raw', 'NA', photoperiod_h, bandNames, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTS)]; %#ok<AGROW>

                    ridgePhaseRows = [ridgePhaseRows; append_ridge_phase_rows( ...
                        job.fileStem, signalID, 'Raw', 'NA', photoperiod_h, bandNames, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTS)]; %#ok<AGROW>

                    bandColourRGB = get_band_colours_rgb(bandNames, BAND_COLOUR);
                    plot_bandridge_timeseries_figs(job.fileStem, signalID, 'Raw', time_day, condChangeIdx, bandNames, bandTS, ...
                        dirs, NAMES, SAVE_DPI, FORCE_ZERO_Y_BP_RPOW, FORCE_ZERO_Y_RP, bandColourRGB);

                    % wide
                    for b = 1:numel(bandNames)
                        wideCol = matlab.lang.makeValidName([signalID '_' bandNames{b}]);
                        wideBP_Raw.(wideCol)   = bandTS.BandPower(b,:).';
                        wideRP_Raw.(wideCol)   = bandTS.RidgePeriod(b,:).';
                        wideRPOW_Raw.(wideCol) = bandTS.RidgePower(b,:).';
                    end
                end
            end
        end

        % ---------------- HSUB ----------------
        if HS.enabled && HS.mapOK && doHSub && rec.found
            residPath = resolve_residual_workbook_path(rec.ResidualWorkbookRel, HS.residualRoot, HS.residualFiles);

            if isempty(residPath) || ~isfile(residPath)
                log_line(RUNLOG, '  HSub skipped: residual workbook not found for "%s"', rec.ResidualWorkbookRel);
            else
                [xH, okHS, errHS] = load_residual_signal_from_workbook(residPath, rec.ResidualColumn, residCache);
                if ~okHS
                    log_line(RUNLOG, '  HSub skipped: %s', errHS);
                else
                    xH = standardise_signal_to_length(xH, N);
                    xH(~isfinite(xH)) = 0;

                    [specOutH, peaksOutH, plotMetaAdd, detailCounter, detailIndexRows] = process_one_stream( ...
                        xH, FB, N, time_day, condChangeIdx, ...
                        job.fileStem, signalID, 'HSub', dirs, NAMES, SAVE_DPI, ...
                        lightStateNumVec, photoperiod_h, ...
                        SPEC_TEMPLATE, SCALO_YTICKS_H, SCALO_HLINES_H, ...
                        BAND_CR, BANDS_UR, detailCacheDir, detailCounter, detailIndexRows, ...
                        BAND_COLOUR);

                    if ~isempty(plotMetaAdd), plotMetaRows = [plotMetaRows; plotMetaAdd]; end %#ok<AGROW>

                    streamStoreHSub(c) = specOutH;
                    if specOutH.ok
                        globalPowerMax_HSub = max(globalPowerMax_HSub, specOutH.localPowerMax);
                    end
                    if ~isempty(peaksOutH), peakRows = [peakRows; peaksOutH]; end %#ok<AGROW>

                    if specOutH.ok
                        bandRowH = compute_band_summary_row(job.fileStem, signalID, 'HSub', photoperiod_h, specOutH.periods_hours, specOutH.avgPowerSpectrum, BAND_CR, BANDS_UR);
                        if ~isempty(bandRowH), bandRows(end+1,:) = bandRowH; end %#ok<SAGROW>

                        bandNamesH = BAND_NAMES_UR;
                        if isfield(specOutH,'defaultBandTS') && ~isempty(specOutH.defaultBandTS)
                            bandTSH = specOutH.defaultBandTS;
                        else
                            bandTSH = compute_band_timeseries(specOutH.periods_hours, specOutH.logPow, BANDS_UR, specOutH.coi_hours);
                        end
                        hsubModeTag = rec.ResidualMode;
                        if isempty(hsubModeTag), hsubModeTag = HS_PRIMARY_VALIDATION_RESIDUAL; end

                        [bandPowRows, ridgePRows, ridgePowRows] = append_band_timeseries_rows_with_valid( ...
                            bandPowRows, ridgePRows, ridgePowRows, ...
                            job.fileStem, signalID, 'HSub', bandNamesH, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTSH);

                        bandCondRows = [bandCondRows; append_band_condition_summary_linearfrac( ...
                            job.fileStem, signalID, 'HSub', photoperiod_h, phaseMasks, bandNamesH, bandTSH)]; %#ok<AGROW>

                        periodCandidateRows = [periodCandidateRows; append_period_candidate_rows( ...
                            job.fileStem, signalID, 'HSub', hsubModeTag, photoperiod_h, phaseMasks, bandNamesH, bandTSH, ...
                            CANDIDATE_MIN_RIDGE_COVERAGE, CANDIDATE_MIN_COI_VALID_FRAC)]; %#ok<AGROW>

                        ridgeTrajectoryRows = [ridgeTrajectoryRows; append_ridge_trajectory_rows( ...
                            job.fileStem, signalID, 'HSub', hsubModeTag, photoperiod_h, bandNamesH, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTSH)]; %#ok<AGROW>

                        ridgePhaseRows = [ridgePhaseRows; append_ridge_phase_rows( ...
                            job.fileStem, signalID, 'HSub', hsubModeTag, photoperiod_h, bandNamesH, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTSH)]; %#ok<AGROW>

                        bandColourRGBH = get_band_colours_rgb(bandNamesH, BAND_COLOUR);
                        plot_bandridge_timeseries_figs(job.fileStem, signalID, 'HSub', time_day, condChangeIdx, bandNamesH, bandTSH, ...
                            dirs, NAMES, SAVE_DPI, FORCE_ZERO_Y_BP_RPOW, FORCE_ZERO_Y_RP, bandColourRGBH);

                        for b = 1:numel(bandNamesH)
                            wideCol = matlab.lang.makeValidName([signalID '_' bandNamesH{b}]);
                            wideBP_HSub.(wideCol)   = bandTSH.BandPower(b,:).';
                            wideRP_HSub.(wideCol)   = bandTSH.RidgePeriod(b,:).';
                            wideRPOW_HSub.(wideCol) = bandTSH.RidgePower(b,:).';
                        end
                    end
                end
            end
        end
    end

    globalPowerMax_Raw  = finalise_global_power_max(globalPowerMax_Raw);
    globalPowerMax_HSub = finalise_global_power_max(globalPowerMax_HSub);

    %% ----------------------- SECOND PASS: POWER PLOTS --------------------
    for c = 1:numel(dataIdx)
        if doRaw
            S = streamStoreRaw(c);
            if S.ok
                plot_power_plots_for_stream(S, job.fileStem, 'Raw', dirs, NAMES, SAVE_DPI, globalPowerMax_Raw, GLOBAL_YMAX_HR);
            end
        end
        if HS.enabled && doHSub
            S = streamStoreHSub(c);
            if S.ok
                plot_power_plots_for_stream(S, job.fileStem, 'HSub', dirs, NAMES, SAVE_DPI, globalPowerMax_HSub, GLOBAL_YMAX_HR);
            end
        end
    end

    %% ------------------- Per-file averaged Abs scalograms ----------------
    % REQUIRED:
    %  - JPEG only (no MAT files)
    %  - identical CLim scaling between Raw and HSub for the same group (p99 across both)
    try
        avgRaw = struct(); avgRaw.ok = false;
        avgH   = struct(); avgH.ok = false;

        if doRaw
            avgRaw = compute_avg_abs_scalograms_perfile(streamStoreRaw, condMap, job.fileStem, 'Raw');
        end
        if HS.enabled && doHSub
            avgH = compute_avg_abs_scalograms_perfile(streamStoreHSub, condMap, job.fileStem, 'HSub');
        end

        % Save per-file averaged plots (Raw/HSub scaling coupled by group where possible)
        save_perfile_avgabs_jpegs(avgRaw, avgH, dirs, NAMES, SAVE_DPI);

        % Add to cross-file pool (store sums/counts only for later combined JPEGs)
        if strcmpi(runMode,'Multiple files') && DO_CROSSFILE_POOLED_AVG
            if doRaw && isfield(avgRaw,'ok') && avgRaw.ok
                CrossPool = crosspool_add_file(CrossPool, avgRaw, photoperiod_h, 'Raw');
            end
            if HS.enabled && doHSub && isfield(avgH,'ok') && avgH.ok
                CrossPool = crosspool_add_file(CrossPool, avgH, photoperiod_h, 'HSub');
            end
        end
    catch ME
        log_line(RUNLOG, 'Averaged Abs scalograms failed: %s', ME.message);
    end

    %% -------------------------- BUILD TABLES -----------------------------
    InputsTable = cell2table(inputsRows, 'VariableNames', inputsHdr);

    PeaksTable = table();
    if ~isempty(peakRows), PeaksTable = cell2table(peakRows, 'VariableNames', peakHdr); end

    DetailIndexTable = table();
    if ~isempty(detailIndexRows), DetailIndexTable = cell2table(detailIndexRows, 'VariableNames', detailIndexHdr); end

    BandSummaryTable = table();
    if ~isempty(bandRows), BandSummaryTable = cell2table(bandRows, 'VariableNames', bandHdr); end

    BandConditionSummaryTable = table();
    if ~isempty(bandCondRows), BandConditionSummaryTable = cell2table(bandCondRows, 'VariableNames', bandCondHdr); end

    PlotMetaTable = table();
    if ~isempty(plotMetaRows), PlotMetaTable = cell2table(plotMetaRows, 'VariableNames', plotMetaHdr); end

    PeriodCandidatesTable = rows_to_table_or_empty(periodCandidateRows, periodCandidateHdr);
    RidgeTrajectoryTable  = rows_to_table_or_empty(ridgeTrajectoryRows,  ridgeTrajectoryHdr);
    RidgePhaseTable       = rows_to_table_or_empty(ridgePhaseRows,       ridgePhaseHdr);

    %% ----------------------------- WRITE EXCEL ---------------------------
    write_summary_workbook(summaryXLSX, InputsTable, PeaksTable, BandSummaryTable, BandConditionSummaryTable, ...
        PeriodCandidatesTable, PERIOD_MIN_MIN, PERIOD_MAX_MIN, TOPN, runType, DO_PHASE_SPLIT, ...
        HS_PRIMARY_VALIDATION_RESIDUAL, CANDIDATE_MIN_RIDGE_COVERAGE, CANDIDATE_MIN_COI_VALID_FRAC);

    write_detail_workbook(detailXLSX, DetailIndexTable, detailCacheDir);

    write_band_timeseries_workbook(bandTSXLSX, bandPowRows, ridgePRows, ridgePowRows, bandTSHdr, RidgeTrajectoryTable, RidgePhaseTable);

    if ~isempty(PlotMetaTable)
        writetable(PlotMetaTable, plotMetaXLSX, 'Sheet', 'PlotMeta');
    end

    %% --------------------- WRITE WIDE TIME SERIES ------------------------
    wideOut = fullfile(dirs.WideTS, sprintf('WideTS_%s.xlsx', sanitise_filename(job.fileStem)));
    delete_if_exists(wideOut);

    writecell({'Wide time-series export per file.'}, wideOut, 'Sheet', 'README');
    writetable(wideBP_Raw,    wideOut, 'Sheet', 'BandPower_Raw');
    writetable(wideBP_HSub,   wideOut, 'Sheet', 'BandPower_HSub');
    writetable(wideRP_Raw,    wideOut, 'Sheet', 'RidgePeriod_Raw');
    writetable(wideRP_HSub,   wideOut, 'Sheet', 'RidgePeriod_HSub');
    writetable(wideRPOW_Raw,  wideOut, 'Sheet', 'RidgePower_Raw');
    writetable(wideRPOW_HSub, wideOut, 'Sheet', 'RidgePower_HSub');

    %% ----------------------------- HANDOFF: SUMMARY ----------------------
    pkgS = struct();
    pkgS.meta = struct();
    pkgS.meta.Script = 'Behav_wavelet_v12_validated_raw_UR_ridge_phase_handoff';
    pkgS.meta.Timestamp = datestr(now,31);
    pkgS.meta.RunType = runType;
    pkgS.meta.doRaw = doRaw;
    pkgS.meta.doHSub = doHSub;
    pkgS.meta.HSubMode = HS_MODE;
    pkgS.meta.PhaseSplitEnabled = DO_PHASE_SPLIT;
    pkgS.meta.PeriodLimits_min = [PERIOD_MIN_MIN, PERIOD_MAX_MIN];
    pkgS.meta.BAND_CR = BAND_CR;
    pkgS.meta.BANDS_UR = BANDS_UR;
    pkgS.meta.BANDS_RAW_ALL = BANDS_RAW_ALL;
    pkgS.meta.BAND_NAMES_ALL = BAND_NAMES_ALL;
    pkgS.meta.BAND_NAMES_UR = BAND_NAMES_UR;
    pkgS.meta.HS_PRIMARY_VALIDATION_RESIDUAL = HS_PRIMARY_VALIDATION_RESIDUAL;
    pkgS.meta.HS_SECONDARY_VALIDATION_RESIDUAL = HS_SECONDARY_VALIDATION_RESIDUAL;
    pkgS.meta.CANDIDATE_MIN_RIDGE_COVERAGE = CANDIDATE_MIN_RIDGE_COVERAGE;
    pkgS.meta.CANDIDATE_MIN_COI_VALID_FRAC = CANDIDATE_MIN_COI_VALID_FRAC;

    pkgS.file = struct();
    pkgS.file.FileStem = string(job.fileStem);
    pkgS.file.InputFile = string(job.inputFile);
    pkgS.file.PerFileFolder = string(job.outDir);

    pkgS.tables = struct();
    pkgS.tables.Inputs = InputsTable;
    pkgS.tables.Peaks = PeaksTable;
    pkgS.tables.DetailIndex = DetailIndexTable;
    pkgS.tables.BandSummary = BandSummaryTable;
    pkgS.tables.BandConditionSummary = BandConditionSummaryTable;
    pkgS.tables.PlotMeta = PlotMetaTable;
    pkgS.tables.PeriodCandidates_Long = PeriodCandidatesTable;

    summaryMat = fullfile(handoffDir, sprintf(NAMES.SUM_PREFIX, sanitise_filename(job.fileStem)));
    save(summaryMat, 'pkgS', '-v7.3');

    %% ----------------------------- HANDOFF: TS ---------------------------
    pkgTS = struct();
    pkgTS.meta = pkgS.meta;
    pkgTS.file = pkgS.file;
    pkgTS.tables = struct();
    pkgTS.tables.BandPower_Long   = rows_to_table_or_empty(bandPowRows, bandTSHdr);
    pkgTS.tables.RidgePeriod_Long = rows_to_table_or_empty(ridgePRows,  bandTSHdr);
    pkgTS.tables.RidgePower_Long  = rows_to_table_or_empty(ridgePowRows, bandTSHdr);
    pkgTS.tables.RidgeTrajectory_Long = RidgeTrajectoryTable;
    pkgTS.tables.RidgePhase_Long = RidgePhaseTable;

    tsMat = fullfile(handoffDir, sprintf(NAMES.TS_PREFIX, sanitise_filename(job.fileStem)));
    save(tsMat, 'pkgTS', '-v7.3');

    handoffIndexRows(end+1,:) = {job.fileStem, job.inputFile, job.outDir, summaryMat, tsMat, runType, HS_MODE, DO_PHASE_SPLIT}; %#ok<SAGROW>

    log_line(RUNLOG, 'Per-file complete. SummaryMat=%s | TSMat=%s', summaryMat, tsMat);
    log_close(RUNLOG);
end

%% ========================================================================
% 9) Cross-file combined averages (Multiple files only)
% ========================================================================
% REQUIRED:
%  - only run when runMode==Multiple files and DO_CROSSFILE_POOLED_AVG true
%  - write single combined scalogram per group
%  - concatenate photoperiod blocks along time axis
%  - boundary lines: white, dotted, thickness comparable
%  - enforce identical CLim between Raw and HSub for same group (p99 across both)
%  - JPEG only (no MAT)
if strcmpi(runMode,'Multiple files') && DO_CROSSFILE_POOLED_AVG
    try
        write_crossfile_combined_abs_jpegs(crossCombinedRoot, CrossPool, SAVE_DPI);
    catch ME
        fprintf('Cross-file combined averages failed: %s\n', ME.message);
    end
end

%% ========================================================================
% 10) Write index
% ========================================================================
idxXLSX = fullfile(handoffDir, NAMES.HANDOFF_INDEX);
delete_if_exists(idxXLSX);

writecell({'Index of per-file handoff mats for AcrossPhotoperiod script.'}, idxXLSX, 'Sheet', 'README');
if ~isempty(handoffIndexRows)
    Tidx = cell2table(handoffIndexRows, 'VariableNames', handoffIndexHdr);
    writetable(Tidx, idxXLSX, 'Sheet', 'Index');
end

fprintf('\nAll done.\nParent output:\n  %s\nHandoff:\n  %s\n', parentOut, handoffDir);

%% ========================================================================
% Local functions
% ========================================================================

function delete_if_exists(p)
    if exist(p,'file'), delete(p); end
end

function dirs = build_perfile_dirs(perFileOutDir, NAMES)
    dirs = struct();
    dirs.Root     = perFileOutDir;
    dirs.Reports  = fullfile(perFileOutDir, NAMES.REPORTS_DIR);
    dirs.FigRoot  = fullfile(perFileOutDir, NAMES.FIG_DIR);
    dirs.Logs     = fullfile(dirs.Reports, NAMES.LOGS_DIR);
    dirs.WideTS   = fullfile(dirs.Reports, NAMES.WIDE_DIR);
    dirs.DetailCache = fullfile(dirs.Reports, NAMES.DETAILCACHE_DIR);

    % Wavelet
    dirs.FigWaveletRoot = fullfile(dirs.FigRoot, NAMES.FIG_WAVELET);
    dirs.FigWavAbs      = fullfile(dirs.FigWaveletRoot, NAMES.FIG_WAV_ABS);
    dirs.FigWavNorm     = fullfile(dirs.FigWaveletRoot, NAMES.FIG_WAV_NORM);
    dirs.FigWavBands    = fullfile(dirs.FigWaveletRoot, NAMES.FIG_WAV_BANDS);
    dirs.FigWavRidge    = fullfile(dirs.FigWaveletRoot, NAMES.FIG_WAV_RIDGE);

    % Averages
    dirs.FigWavAvgRoot  = fullfile(dirs.FigWaveletRoot, NAMES.FIG_WAV_AVG);
    dirs.FigWavAvgRaw   = fullfile(dirs.FigWavAvgRoot, 'Raw');
    dirs.FigWavAvgHSub  = fullfile(dirs.FigWavAvgRoot, 'HSub');

    % Spectra
    dirs.FigPowerRoot = fullfile(dirs.FigRoot, NAMES.FIG_POWER);
    dirs.FigPwr_Raw_Glob  = fullfile(dirs.FigPowerRoot, NAMES.PWR_RAW,  NAMES.PWR_GLOBAL);
    dirs.FigPwr_Raw_Cond  = fullfile(dirs.FigPowerRoot, NAMES.PWR_RAW,  NAMES.PWR_COND);
    dirs.FigPwr_H_Glob    = fullfile(dirs.FigPowerRoot, NAMES.PWR_HSUB, NAMES.PWR_GLOBAL);
    dirs.FigPwr_H_Cond    = fullfile(dirs.FigPowerRoot, NAMES.PWR_HSUB, NAMES.PWR_COND);

    % Band/Ridge
    dirs.FigBandRidgeRoot = fullfile(dirs.FigRoot, NAMES.FIG_BANDRIDGE);
    dirs.BR_BP_Comb    = fullfile(dirs.FigBandRidgeRoot, NAMES.BR_BP,   NAMES.BR_COMBINED);
    dirs.BR_RP_Comb    = fullfile(dirs.FigBandRidgeRoot, NAMES.BR_RP,   NAMES.BR_COMBINED);
    dirs.BR_RPOW_Comb  = fullfile(dirs.FigBandRidgeRoot, NAMES.BR_RPOW, NAMES.BR_COMBINED);
end

function ensure_perfile_dirs(dirs)
    ensure_dir(dirs.Reports);
    ensure_dir(dirs.Logs);
    ensure_dir(dirs.FigRoot);
    ensure_dir(dirs.WideTS);
    ensure_dir(dirs.DetailCache);

    ensure_dir(dirs.FigWaveletRoot);
    ensure_dir(dirs.FigWavAbs);
    ensure_dir(dirs.FigWavNorm);
    ensure_dir(dirs.FigWavBands);
    ensure_dir(dirs.FigWavRidge);

    ensure_dir(dirs.FigWavAvgRoot);
    ensure_dir(dirs.FigWavAvgRaw);
    ensure_dir(dirs.FigWavAvgHSub);

    ensure_dir(dirs.FigPowerRoot);
    ensure_dir(dirs.FigPwr_Raw_Glob);
    ensure_dir(dirs.FigPwr_Raw_Cond);
    ensure_dir(dirs.FigPwr_H_Glob);
    ensure_dir(dirs.FigPwr_H_Cond);

    ensure_dir(dirs.FigBandRidgeRoot);
    ensure_dir(dirs.BR_BP_Comb);
    ensure_dir(dirs.BR_RP_Comb);
    ensure_dir(dirs.BR_RPOW_Comb);
end

function T = rows_to_table_or_empty(rows, hdr)
    if isempty(rows)
        T = table();
    else
        T = cell2table(rows, 'VariableNames', hdr);
    end
end

function [isLight, lightStateStr, lightStateNum] = derive_light_state_from_lightdur(time_hr, lightDurVec)
    t = double(time_hr(:));
    L = double(lightDurVec(:));
    n = min(numel(t), numel(L));
    t = t(1:n);
    L = L(1:n);

    ZT = mod(t,24);
    isLight = false(n,1);
    ok = isfinite(ZT) & isfinite(L) & L>=0 & L<=24;
    isLight(ok) = (ZT(ok) < L(ok)); % boundary at ZT=L is Dark

    lightStateNum = double(isLight); % Light=1, Dark=0
    ls = repmat("Dark", n, 1);
    ls(isLight) = "Light";
    lightStateStr = cellstr(ls);
end

function [fileList, inputParent] = select_input_files(runMode)
    fileList = {};
    inputParent = '';

    if strcmpi(runMode, 'Single file')
        [fileName, filePath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, 'Select input behavioural Excel file');
        if isequal(fileName,0), return; end
        fileList = {fullfile(filePath, fileName)};
        inputParent = filePath;
    else
        [files, filePath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, 'Select input behavioural Excel file(s)', 'MultiSelect','on');
        if isequal(files,0), return; end
        if ischar(files), files = {files}; end
        files = files(:);
        fileList = cell(numel(files),1);
        for i = 1:numel(files)
            fileList{i} = fullfile(filePath, files{i});
        end
        inputParent = filePath;
    end
end

function ensure_dir(p)
    if ~exist(p, 'dir'), mkdir(p); end
end

function s = sanitise_filename(strIn)
    s = regexprep(string(strIn), '[^\w\-]', '_');
    s = char(s);
    if isempty(s), s = 'x'; end
    if numel(s) > 100, s = s(1:100); end
end

function L = log_open(path)
    L = struct();
    L.path = path;
    L.fid = fopen(path, 'a');
    if L.fid < 0
        warning('Could not open log file: %s', path);
        L.fid = [];
    end
end

function log_line(L, fmt, varargin)
    try
        if isempty(L) || ~isfield(L,'fid') || isempty(L.fid) || L.fid < 0
            return;
        end
        msg = sprintf(fmt, varargin{:});
        fprintf(L.fid, '[%s] %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg);
    catch
    end
end

function log_close(L)
    try
        if ~isempty(L) && isfield(L,'fid') && ~isempty(L.fid) && L.fid > 0
            fclose(L.fid);
        end
    catch
    end
end

%% ----------------------------- HS INDEX ---------------------------------
function HSIndex = build_hsub_index_from_parent(parentFolder)
    HSIndex = struct();
    HSIndex.ok = false;
    HSIndex.byStem = struct();

    if ~exist(parentFolder,'dir'), return; end

    dd = dir(fullfile(parentFolder, '**', 'Reports', 'HS_Summary.xlsx'));
    if isempty(dd), return; end

    cand = struct('stem',{},'hsSummaryPath',{},'runFolder',{},'residualFiles',{},'mtime',{});

    for i = 1:numel(dd)
        hsPath = fullfile(dd(i).folder, dd(i).name);
        stem = parse_hsub_summary_stem_only(hsPath);
        if isempty(stem), continue; end

        runFolder = fileparts(dd(i).folder);
        resFolder = fullfile(runFolder, 'TS', 'RES');
        resFiles = {};
        if exist(resFolder,'dir')
            rf = dir(fullfile(resFolder, '*.xlsx'));
            resFiles = cell(numel(rf),1);
            for k = 1:numel(rf)
                resFiles{k} = fullfile(rf(k).folder, rf(k).name);
            end
        end

        c = struct();
        c.stem = stem;
        c.hsSummaryPath = hsPath;
        c.runFolder = runFolder;
        c.residualFiles = resFiles;
        c.mtime = get_file_mtime(hsPath);
        cand(end+1) = c; %#ok<AGROW>
    end

    if isempty(cand), return; end

    stems = string({cand.stem});
    u = unique(stems);

    for i = 1:numel(u)
        st = char(u(i));
        idx = find(stems == u(i));
        group = cand(idx);
        mt = [group.mtime];
        [~, ord] = sort(mt, 'descend');
        chosen = group(ord(1));

        key = matlab.lang.makeValidName(st);
        HSIndex.byStem.(key) = chosen;
    end

    HSIndex.ok = true;
end

function t = get_file_mtime(p)
    try
        d = dir(p);
        if isempty(d), t = 0; else, t = d(1).datenum; end %#ok<DATNM>
    catch
        t = 0;
    end
end

function stem = parse_hsub_summary_stem_only(hsSummaryPath)
    stem = '';
    if ~isfile(hsSummaryPath), return; end
    try
        sh = sheetnames(hsSummaryPath);
        if isempty(sh), return; end
        raw = readcell(hsSummaryPath, 'Sheet', sh{1});
        if isempty(raw), return; end
        s = string(raw); s = s(:); s(ismissing(s)) = "";
        ix = find(contains(lower(s), "input file"), 1, 'first');
        if isempty(ix), return; end
        line = char(s(ix));
        tok = regexp(line, '([A-Za-z0-9_\-]+)\.(xlsx|xls)', 'match', 'once');
        if ~isempty(tok)
            [~, st] = fileparts(tok);
            stem = st;
        end
    catch
        stem = '';
    end
end

function [entry, found] = hsub_index_lookup(HSIndex, stem)
    entry = struct();
    found = false;
    if isempty(HSIndex) || ~isfield(HSIndex,'ok') || ~HSIndex.ok, return; end
    key = matlab.lang.makeValidName(char(string(stem)));
    if isfield(HSIndex.byStem, key)
        entry = HSIndex.byStem.(key);
        found = true;
    end
end

%% ----------------------------- IO ---------------------------------------
function [dataTable, info] = read_input_table_preserve_robust(inputFile)
    info = struct(); info.File = inputFile; info.Method = '';
    try
        opts = detectImportOptions(inputFile, 'FileType', 'spreadsheet');
        if isprop(opts, 'VariableNamingRule'), opts.VariableNamingRule = 'preserve'; end
        try
            opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA','N/A','',' '});
        catch
        end
        dataTable = readtable(inputFile, opts);
        info.Method = 'detectImportOptions+readtable';
        return;
    catch
        dataTable = readtable(inputFile, 'VariableNamingRule', 'preserve');
        info.Method = 'readtable_preserve_fallback';
    end
end

function [T, info] = drop_empty_columns_robust(T)
    info = struct(); info.DroppedIdx = []; info.DroppedNames = {};
    if isempty(T) || width(T) == 0, return; end
    vn = T.Properties.VariableNames;
    dropMask = false(1, width(T));
    for i = 1:width(T)
        x = T{:, i};
        dropMask(i) = is_column_all_missing(x);
    end
    if any(dropMask)
        info.DroppedIdx = find(dropMask);
        info.DroppedNames = vn(dropMask);
        T(:, dropMask) = [];
    end
end

function tfCol = is_column_all_missing(x)
    try
        if isdatetime(x), tfCol = all(isnat(x)); return; end
        if isnumeric(x) || islogical(x), tfCol = all(isnan(double(x))); return; end
        if isstring(x), tfCol = all(ismissing(x)); return; end
        if iscell(x)
            tfCol = true;
            for k = 1:numel(x)
                v = x{k};
                if isempty(v), continue; end
                if (ischar(v) || (isstring(v)&&isscalar(v))) && strlength(string(v))==0, continue; end
                if isnumeric(v) && isscalar(v) && isnan(v), continue; end
                tfCol = false; return;
            end
            return;
        end
        tfCol = all(ismissing(string(x)));
    catch
        tfCol = false;
    end
end

%% ----------------------------- PREFLIGHT UI ------------------------------
function mapPF = preflight_column_mapping_dialog_time_lightduration_and_conditions(dataTable, fileIndex, nFiles)
    varNames = dataTable.Properties.VariableNames;
    nVars = numel(varNames);

    [listStr, exactKeys] = make_list_labels(varNames);
    titleBar = sprintf('Column setup | File %d/%d', fileIndex, nFiles);

    timeGuess = guess_time_column(varNames);
    if isempty(timeGuess), timeGuess = 1; end
    timeGuess = max(1, min(nVars, timeGuess));

    [timeIdx, okT] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', {sprintf('File %d/%d', fileIndex, nFiles), 'Select Time column'}, ...
        'SelectionMode', 'single', ...
        'ListString', listStr, ...
        'InitialValue', timeGuess, ...
        'ListSize', [720 420]);
    if ~okT || isempty(timeIdx)
        error('No Time column selected. Preflight cancelled.');
    end

    ldGuess = guess_lightduration_column(varNames);
    if isempty(ldGuess), ldGuess = 1; end
    ldGuess = max(1, min(nVars, ldGuess));

    [ldIdx, okLD] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', {sprintf('File %d/%d', fileIndex, nFiles), 'Select Light duration (h) column (required)'}, ...
        'SelectionMode', 'single', ...
        'ListString', listStr, ...
        'InitialValue', ldGuess, ...
        'ListSize', [720 420]);
    if ~okLD || isempty(ldIdx)
        error('Light duration selection cancelled.');
    end

    defaultData = setdiff(1:nVars, [timeIdx, ldIdx], 'stable');
    [dataIdx, okD] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', {sprintf('File %d/%d', fileIndex, nFiles), 'Select Data columns'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listStr, ...
        'InitialValue', defaultData, ...
        'ListSize', [720 420]);
    if ~okD || isempty(dataIdx)
        error('No Data columns selected. Preflight cancelled.');
    end

    dataIdx = setdiff(dataIdx, [timeIdx, ldIdx], 'stable');
    if isempty(dataIdx)
        error('After excluding Time/Light duration, no Data columns remain.');
    end

    selectedDataNames = varNames(dataIdx);
    selectedDataLabels = listStr(dataIdx);
    excludeList = [{'<<None>>'}; selectedDataLabels(:)];
    [excludeSel, okE] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', {sprintf('File %d/%d', fileIndex, nFiles), 'Optional: select Data columns to exclude (or keep <<None>>)'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', excludeList, ...
        'InitialValue', 1, ...
        'ListSize', [720 420]);

    excludedNames = {};
    if okE
        excludeSel = setdiff(excludeSel, 1);
        if ~isempty(excludeSel)
            exPos = excludeSel - 1;
            exPos = exPos(exPos >= 1 & exPos <= numel(dataIdx));
            if ~isempty(exPos)
                excludedNames = selectedDataNames(exPos);
                dataIdx(exPos) = [];
            end
        end
    end

    if isempty(dataIdx)
        error('All Data columns were excluded.');
    end

    % Conditions-by-columns UI
    dataNames = varNames(dataIdx);
    nData = numel(dataNames);

    ansN = inputdlg({sprintf('Number of conditions for this file (1..12):')}, ...
        sprintf('Conditions | File %d/%d', fileIndex, nFiles), [1 60], {'2'});
    if isempty(ansN), error('Condition setup cancelled.'); end
    nCond = round(str2double(ansN{1}));
    if ~isfinite(nCond) || nCond < 1 || nCond > 12
        error('Invalid number of conditions.');
    end

    condLabels = cell(nCond,1);
    for j = 1:nCond
        a = inputdlg({sprintf('Label for Condition %d:', j)}, ...
            sprintf('Conditions | File %d/%d', fileIndex, nFiles), [1 60], {sprintf('Cond%d',j)});
        if isempty(a), error('Condition label entry cancelled.'); end
        lab = char(strtrim(string(a{1})));
        if isempty(lab), lab = sprintf('Cond%d',j); end
        condLabels{j} = lab;
    end

    remaining = dataNames(:);
    assigned = containers.Map('KeyType','char','ValueType','char');

    for j = 1:nCond
        if isempty(remaining), break; end
        listD = cellstr(string(remaining));
        [sel, ok] = listdlg( ...
            'Name', sprintf('Assign columns | %s', condLabels{j}), ...
            'PromptString', {sprintf('Assign columns to "%s"', condLabels{j}), 'Select one or more columns'}, ...
            'SelectionMode','multiple', ...
            'ListString', listD, ...
            'ListSize', [720 420]);
        if ~ok
            error('Condition assignment cancelled.');
        end

        chosen = remaining(sel);
        for k = 1:numel(chosen)
            assigned(char(string(chosen(k)))) = condLabels{j};
        end

        keepMask = true(numel(remaining),1);
        keepMask(sel) = false;
        remaining = remaining(keepMask);
    end

    if ~isempty(remaining)
        for k = 1:numel(remaining)
            assigned(char(string(remaining(k)))) = 'Unassigned';
        end
    end

    mapPF = struct();
    mapPF.nColsPreflight = nVars;
    mapPF.varNamesPreflight = varNames(:)';
    mapPF.listLabelsPreflight = listStr(:)';
    mapPF.listExactKeys = exactKeys(:)';

    mapPF.timeIdx_pre = timeIdx;
    mapPF.lightDurIdx_pre = ldIdx;

    mapPF.timeName_pre = varNames{timeIdx};
    mapPF.timeCanon = canonicalise_header(varNames{timeIdx});

    mapPF.lightDurName_pre = varNames{ldIdx};
    mapPF.lightDurCanon = canonicalise_header(varNames{ldIdx});

    mapPF.dataNames_pre = varNames(dataIdx);
    mapPF.dataCanon = cellfun(@canonicalise_header, varNames(dataIdx), 'UniformOutput', false);
    mapPF.excludedNames = excludedNames(:)';

    mapPF.condLabels = condLabels(:)';
    mapPF.condAssigned = assigned;
    mapPF.nCond = nCond;
    mapPF.nData = nData;
end

function [labels, exactKeys] = make_list_labels(varNames)
    n = numel(varNames);
    labels = cell(n,1);
    exactKeys = cell(n,1);

    raw = cell(n,1);
    for i = 1:n
        raw{i} = char(string(varNames{i}));
        if isempty(strtrim(raw{i}))
            raw{i} = sprintf('(blank_header_col_%d)', i);
        end
    end

    [u, ~, ic] = unique(raw, 'stable');
    counts = accumarray(ic, 1);
    dupCounter = zeros(size(u));

    for i = 1:n
        label = raw{i};
        g = ic(i);
        if counts(g) > 1
            dupCounter(g) = dupCounter(g) + 1;
            label = sprintf('%s  [#%d]', label, dupCounter(g));
        end
        labels{i} = label;
        exactKeys{i} = raw{i};
    end
end

function idx = guess_time_column(varNames)
    idx = [];
    canons = cellfun(@canonicalise_header, varNames, 'UniformOutput', false);
    strong = find(cellfun(@(s) contains(s,'time') && (contains(s,'hr') || contains(s,'hour') || contains(s,'min') || contains(s,'date')), canons), 1, 'first');
    if ~isempty(strong), idx = strong; return; end
    anyTime = find(cellfun(@(s) contains(s,'time') || contains(s,'datetime') || contains(s,'date'), canons), 1, 'first');
    if ~isempty(anyTime), idx = anyTime; end
end

function idx = guess_lightduration_column(varNames)
    idx = [];
    canons = cellfun(@canonicalise_header, varNames, 'UniformOutput', false);
    m = find(cellfun(@(s) contains(s,'light') && contains(s,'duration'), canons), 1, 'first');
    if ~isempty(m), idx = m; return; end
    m = find(cellfun(@(s) contains(s,'photoperiod') || contains(s,'lightdur'), canons), 1, 'first');
    if ~isempty(m), idx = m; end
end

function c = canonicalise_header(x)
    c = lower(char(string(x)));
    c = regexprep(c, '\s+', '');
    c = regexprep(c, '[^a-z0-9]', '');
end

function [timeIdx, lightDurIdx, dataIdx, excludedNames, notes, condMap] = resolve_preflight_mapping_time_lightdur_and_conditions(dataTable, mapPF)
    varNamesNow = dataTable.Properties.VariableNames;
    notes = {};

    [timeIdx, noteT] = resolve_one_column(varNamesNow, mapPF.timeName_pre, mapPF.timeCanon, mapPF.timeIdx_pre, 'Time');
    if ~isempty(noteT), notes{end+1} = noteT; end %#ok<AGROW>

    [ldIdx, noteLD] = resolve_one_column(varNamesNow, mapPF.lightDurName_pre, mapPF.lightDurCanon, mapPF.lightDurIdx_pre, 'Light duration (h)');
    lightDurIdx = ldIdx;
    if ~isempty(noteLD), notes{end+1} = noteLD; end %#ok<AGROW>

    dataIdx = [];
    missingData = {};
    nData = numel(mapPF.dataNames_pre);
    canonNow = cellfun(@canonicalise_header, varNamesNow, 'UniformOutput', false);

    for i = 1:nData
        thisName = mapPF.dataNames_pre{i};
        thisCanon = mapPF.dataCanon{i};

        exactMatch = find(strcmp(varNamesNow, thisName), 1, 'first');
        if ~isempty(exactMatch)
            dataIdx(end+1) = exactMatch; %#ok<AGROW>
            continue;
        end

        canonMatch = find(strcmp(canonNow, thisCanon), 1, 'first');
        if ~isempty(canonMatch)
            dataIdx(end+1) = canonMatch; %#ok<AGROW>
            notes{end+1} = sprintf('Data[%d] resolved by canonical match: "%s" -> "%s"', i, thisName, varNamesNow{canonMatch}); %#ok<AGROW>
            continue;
        end

        missingData{end+1} = thisName; %#ok<AGROW>
    end

    dataIdx = unique(dataIdx, 'stable');
    dataIdx = setdiff(dataIdx, [timeIdx, lightDurIdx], 'stable');

    if ~isempty(missingData)
        notes{end+1} = sprintf('Preflight Data columns not found (skipped): %s', strjoin(missingData, ', '));
    end

    if isempty(timeIdx), error('Preflight Time column not found.'); end
    if isempty(lightDurIdx), error('Preflight Light duration column not found.'); end
    if isempty(dataIdx), error('No preflight Data columns could be resolved in the current file.'); end

    excludedNames = {};
    if isfield(mapPF,'excludedNames') && ~isempty(mapPF.excludedNames)
        excludedNames = mapPF.excludedNames;
    end

    % ConditionMap for CURRENT resolved data columns
    dataNamesNow = varNamesNow(dataIdx);
    assigned = containers.Map('KeyType','char','ValueType','char');
    if isfield(mapPF,'condAssigned') && isa(mapPF.condAssigned,'containers.Map')
        assigned = mapPF.condAssigned;
    end

    condOf = cell(numel(dataNamesNow),1);
    for i = 1:numel(dataNamesNow)
        nm = char(string(dataNamesNow{i}));
        if isKey(assigned, nm)
            condOf{i} = assigned(nm);
        else
            condOf{i} = 'Unassigned';
        end
    end
    condMap = table(string(dataNamesNow(:)), string(condOf(:)), 'VariableNames', {'SignalID','ConditionLabel'});
end

function [idx, note] = resolve_one_column(varNamesNow, targetName_pre, targetCanon, idxPre, roleLabel)
    idx = [];
    note = '';

    exactMatch = find(strcmp(varNamesNow, targetName_pre), 1, 'first');
    if ~isempty(exactMatch)
        idx = exactMatch;
        return;
    end

    canonNow = cellfun(@canonicalise_header, varNamesNow, 'UniformOutput', false);
    canonMatch = find(strcmp(canonNow, targetCanon), 1, 'first');
    if ~isempty(canonMatch)
        idx = canonMatch;
        note = sprintf('%s resolved by canonical match: "%s" -> "%s"', roleLabel, targetName_pre, varNamesNow{canonMatch});
        return;
    end

    if ~isempty(idxPre) && isfinite(idxPre) && idxPre >= 1 && idxPre <= numel(varNamesNow)
        idx = idxPre;
        note = sprintf('%s resolved by stored index (%d): "%s"', roleLabel, idxPre, varNamesNow{idxPre});
        return;
    end
end

function write_condition_map_workbook(outXLSX, condMap, fileStem)
    writecell({sprintf('ConditionMap for file: %s', char(string(fileStem))); ...
               'SignalID -> ConditionLabel'; ...
               'Defined during preflight.'}, outXLSX, 'Sheet', 'README');
    if ~isempty(condMap) && istable(condMap)
        writetable(condMap, outXLSX, 'Sheet', 'ConditionMap');
    else
        writecell({'No ConditionMap produced.'}, outXLSX, 'Sheet', 'ConditionMap');
    end
end

%% ----------------------------- TIME -------------------------------------
function [timeMinutes, TsMinutes] = infer_time_minutes(timeCol, timeName)
    if isdatetime(timeCol)
        t = timeCol(:);
        dt = minutes(diff(t));
        TsMinutes = median(dt, 'omitnan');
        timeMinutes = minutes(t - t(1));
        timeMinutes = double(timeMinutes);
        return;
    end

    if isnumeric(timeCol)
        t = double(timeCol(:));
        nm = lower(string(timeName));

        if contains(nm, "min")
            timeMinutes = t;
        elseif contains(nm, "hr") || contains(nm, "hour")
            timeMinutes = t * 60;
        elseif contains(nm, "day")
            timeMinutes = t * 24 * 60;
        else
            dtRaw = diff(t);
            medStep = median(dtRaw, 'omitnan');
            if max(t, [], 'omitnan') > 24 && isfinite(medStep) && medStep > 0 && medStep < 1
                timeMinutes = t * 60;
            else
                timeMinutes = t;
            end
        end

        dt = diff(timeMinutes);
        TsMinutes = median(dt, 'omitnan');
        if ~isfinite(TsMinutes) || TsMinutes <= 0
            error('Sampling interval could not be inferred from the Time column.');
        end

        timeMinutes = timeMinutes - timeMinutes(1);
        return;
    end

    t = str2double(string(timeCol(:)));
    if all(~isfinite(t))
        error('Unsupported Time column type. Time must be numeric or datetime.');
    end
    [timeMinutes, TsMinutes] = infer_time_minutes(t, timeName);
end

function y = to_numeric_vector(x)
    if isnumeric(x) || islogical(x)
        y = double(x(:)); return;
    end
    if isdatetime(x) || isduration(x)
        y = []; return;
    end
    if isstring(x) || ischar(x)
        y = str2double(string(x)); y = y(:); return;
    end
    if iscell(x)
        try
            y = str2double(string(x)); y = y(:);
        catch
            y = [];
        end
        return;
    end
    try
        y = str2double(string(x)); y = y(:);
    catch
        y = [];
    end
end

function x = standardise_signal_to_length(x, N)
    x = to_numeric_vector(x);
    if isempty(x), return; end
    x = x(:);
    if numel(x) < N
        x(end+1:N,1) = NaN; %#ok<AGROW>
    elseif numel(x) > N
        x = x(1:N);
    end
end

function phot = robust_scalar_mode(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x), phot = NaN; return; end
    u = unique(x);
    if numel(u) == 1, phot = u; return; end
    phot = mode(x);
    if ~isfinite(phot), phot = median(x,'omitnan'); end
end

function idxChange = compute_changes_from_lightduration(lightDurVec)
    idxChange = [];
    try
        v = double(lightDurVec(:));
        if isempty(v), return; end
        dv = diff(v);
        idxChange = find(dv ~= 0 & isfinite(dv));
    catch
        idxChange = [];
    end
end

%% ----------------------------- HS MAP -----------------------------------
function [map, ok] = load_hsub_recommend_map(hsSummaryPath)
    % v12 behaviour:
    %   Prefer the ValidationManifest sheet produced by Harmonic_subtraction_v12.
    %   This forces downstream HSub wavelets to use PrimaryValidationResidual
    %   (SEL_P360 by design) while retaining recommendation metadata.
    % Fallback:
    %   If ValidationManifest is absent, use the legacy Recommend sheet.
    map = struct(); ok = false;
    if ~isfile(hsSummaryPath), return; end

    try
        sh = sheetnames(hsSummaryPath);
    catch
        return;
    end

    map.byValidName = struct();
    map.byCanonical = struct();
    map.SourceSheet = '';

    % ---------------- v12 preferred path: ValidationManifest -------------
    if any(strcmpi(sh, 'ValidationManifest'))
        try
            R = readtable(hsSummaryPath, 'Sheet', 'ValidationManifest', 'VariableNamingRule', 'preserve');
        catch
            R = table();
        end

        if ~isempty(R) && height(R) > 0
            vns = string(R.Properties.VariableNames);
            req = {'Column','PrimaryValidationResidual','PrimaryValidationWorkbook','PrimaryValidationSignalColumnHeader'};
            hasReq = true;
            for i = 1:numel(req)
                hasReq = hasReq && any(strcmpi(vns, req{i}));
            end

            if hasReq
                col_Column = get_col_case_insensitive(R, 'Column');
                col_Mode   = get_col_case_insensitive(R, 'PrimaryValidationResidual');
                col_WB     = get_col_case_insensitive(R, 'PrimaryValidationWorkbook');
                col_SigCol = get_col_case_insensitive(R, 'PrimaryValidationSignalColumnHeader');

                col_RecMode = find_col_case_insensitive_optional(R, 'RecommendedResidual');
                col_RecWB   = find_col_case_insensitive_optional(R, 'RecommendedWorkbook');
                col_FullFlg = find_col_case_insensitive_optional(R, 'FullLadderRecommendedFlag');
                col_QCFlag  = find_col_case_insensitive_optional(R, 'HSubQCFlag');

                for r = 1:height(R)
                    rawCol = string(R{r, col_Column});
                    if strlength(rawCol) == 0 || ismissing(rawCol), continue; end
                    rawColChar = char(rawCol);

                    wbRel     = string(R{r, col_WB});
                    sigColHdr = string(R{r, col_SigCol});
                    mode      = string(R{r, col_Mode});

                    entry = struct();
                    entry.ResidualWorkbookRel = char(wbRel);
                    entry.ResidualColumn      = char(sigColHdr);
                    entry.ResidualMode        = char(mode);
                    entry.found               = true;
                    entry.MapSourceSheet      = 'ValidationManifest';
                    entry.RecommendedResidual = '';
                    entry.RecommendedWorkbook = '';
                    entry.FullLadderRecommendedFlag = false;
                    entry.HSubQCFlag = '';

                    if col_RecMode > 0, entry.RecommendedResidual = char(string(R{r, col_RecMode})); end
                    if col_RecWB   > 0, entry.RecommendedWorkbook = char(string(R{r, col_RecWB})); end
                    if col_FullFlg > 0, entry.FullLadderRecommendedFlag = scalar_to_logical(R{r, col_FullFlg}); end
                    if col_QCFlag  > 0, entry.HSubQCFlag = char(string(R{r, col_QCFlag})); end

                    validField = matlab.lang.makeValidName(rawColChar);
                    canonKey   = canonicalise_header(rawColChar);
                    map.byValidName.(validField) = entry;
                    map.byCanonical.(matlab.lang.makeValidName(canonKey)) = entry;
                end
                map.SourceSheet = 'ValidationManifest';
                ok = true;
                return;
            end
        end
    end

    % ---------------- legacy fallback: Recommend --------------------------
    if ~any(strcmpi(sh, 'Recommend')), return; end
    try
        R = readtable(hsSummaryPath, 'Sheet', 'Recommend', 'VariableNamingRule', 'preserve');
    catch
        return;
    end
    if isempty(R) || height(R) == 0, return; end

    req = {'Column','RecommendedResidual','Workbook','SignalColumnHeader'};
    vns = string(R.Properties.VariableNames);
    for i = 1:numel(req)
        if ~any(strcmpi(vns, req{i})), return; end
    end

    col_Column = get_col_case_insensitive(R, 'Column');
    col_Mode   = get_col_case_insensitive(R, 'RecommendedResidual');
    col_WB     = get_col_case_insensitive(R, 'Workbook');
    col_SigCol = get_col_case_insensitive(R, 'SignalColumnHeader');

    for r = 1:height(R)
        rawCol = string(R{r, col_Column});
        if strlength(rawCol) == 0 || ismissing(rawCol), continue; end
        rawColChar = char(rawCol);

        validField = matlab.lang.makeValidName(rawColChar);
        canonKey   = canonicalise_header(rawColChar);

        wbRel    = string(R{r, col_WB});
        sigColHdr= string(R{r, col_SigCol});
        mode     = string(R{r, col_Mode});

        entry = struct();
        entry.ResidualWorkbookRel = char(wbRel);
        entry.ResidualColumn      = char(sigColHdr);
        entry.ResidualMode        = char(mode);
        entry.found               = true;
        entry.MapSourceSheet      = 'Recommend';
        entry.RecommendedResidual = char(mode);
        entry.RecommendedWorkbook = char(wbRel);
        entry.FullLadderRecommendedFlag = startsWith(mode, 'FL');
        entry.HSubQCFlag = '';

        map.byValidName.(validField) = entry;
        map.byCanonical.(matlab.lang.makeValidName(canonKey)) = entry;
    end

    map.SourceSheet = 'Recommend';
    ok = true;
end

function idx = get_col_case_insensitive(T, colName)
    v = string(T.Properties.VariableNames);
    m = find(strcmpi(v, colName), 1, 'first');
    if isempty(m), error('Column not found: %s', colName); end
    idx = m;
end

function idx = find_col_case_insensitive_optional(T, colName)
    v = string(T.Properties.VariableNames);
    m = find(strcmpi(v, colName), 1, 'first');
    if isempty(m), idx = 0; else, idx = m; end
end

function tf = scalar_to_logical(x)
    try
        if iscell(x), x = x{1}; end
        if islogical(x), tf = logical(x); return; end
        if isnumeric(x), tf = isfinite(x) && x ~= 0; return; end
        xs = lower(strtrim(char(string(x))));
        tf = any(strcmp(xs, {'true','1','yes','y'}));
    catch
        tf = false;
    end
end

function rec = get_recommendation_for_column(map, signalID)
    rec = struct(); rec.found = false;
    rec.ResidualWorkbookRel = ''; rec.ResidualColumn = ''; rec.ResidualMode = '';
    rec.MapSourceSheet = ''; rec.RecommendedResidual = ''; rec.RecommendedWorkbook = '';
    rec.FullLadderRecommendedFlag = false; rec.HSubQCFlag = '';

    if isempty(map) || ~isfield(map,'byValidName'), return; end

    fld = matlab.lang.makeValidName(signalID);
    if isfield(map.byValidName, fld)
        rec = map.byValidName.(fld); rec.found = true; return;
    end
    canonKey = canonicalise_header(signalID);
    canonFld = matlab.lang.makeValidName(canonKey);
    if isfield(map.byCanonical, canonFld)
        rec = map.byCanonical.(canonFld); rec.found = true; return;
    end
end

function residPath = resolve_residual_workbook_path(workbookRel, residualRoot, selectedResidualFiles)
    residPath = '';
    rel = char(string(workbookRel));
    if isempty(rel), return; end

    if ~isempty(residualRoot)
        root = char(string(residualRoot));
        rootN = lower(strrep(root,'/','\'));
        rel2 = rel;
        rel2N = lower(strrep(rel2,'/','\'));
        if endsWith(rootN, '\ts\res') && startsWith(rel2N, 'ts\res\')
            rel2 = rel2(8:end);
        elseif endsWith(rootN, '\ts\res\') && startsWith(rel2N, 'ts\res\')
            rel2 = rel2(8:end);
        end
        cand = fullfile(root, rel2);
        if isfile(cand), residPath = cand; return; end
    end

    if ~isempty(selectedResidualFiles)
        relKey = lower(strrep(rel,'/','\'));
        for i = 1:numel(selectedResidualFiles)
            p = lower(strrep(char(string(selectedResidualFiles{i})),'/','\'));
            if endsWith(p, relKey)
                residPath = selectedResidualFiles{i}; return;
            end
        end
    end

    [~, relName, relExt] = fileparts(rel);
    relFile = lower([relName relExt]);
    if ~isempty(selectedResidualFiles)
        for i = 1:numel(selectedResidualFiles)
            [~, nm, ex] = fileparts(selectedResidualFiles{i});
            if strcmpi([nm ex], relFile)
                residPath = selectedResidualFiles{i}; return;
            end
        end
    end
end

function [sig, ok, errMsg] = load_residual_signal_from_workbook(residPath, sigColHdr, residCache)
    sig = []; ok = false; errMsg = '';
    try
        if ~isfile(residPath)
            errMsg = sprintf('Residual workbook not found: %s', residPath);
            return;
        end

        key = lower(char(residPath));
        if isKey(residCache, key)
            residTab = residCache(key);
        else
            residTab = readtable(residPath, 'VariableNamingRule','preserve');
            residCache(key) = residTab;
        end

        vns = residTab.Properties.VariableNames;
        cExact = find(strcmp(vns, sigColHdr), 1, 'first');
        if isempty(cExact)
            canonTarget = canonicalise_header(sigColHdr);
            canonNow = cellfun(@canonicalise_header, vns, 'UniformOutput', false);
            cCanon = find(strcmp(canonNow, canonTarget), 1, 'first');
            if isempty(cCanon)
                errMsg = sprintf('Residual column not found: "%s"', sigColHdr);
                return;
            else
                cExact = cCanon;
            end
        end

        sig = to_numeric_vector(residTab{:, cExact});
        if isempty(sig)
            errMsg = sprintf('Residual signal non-numeric in column "%s".', sigColHdr);
            return;
        end

        ok = true;
    catch ME
        errMsg = sprintf('Residual load failed: %s', ME.message);
    end
end

%% ----------------------------- SPEC TEMPLATE -----------------------------
function S = make_spec_template()
    S = struct();
    S.ok = false;
    S.signalID = '';
    S.sourceTag = '';
    S.periods_hours = [];
    S.time_day = [];
    S.logPow = [];
    S.coi_hours = [];
    S.avgPowerSpectrum = [];
    S.condLabels = {};
    S.condSpectra = [];
    S.localPowerMax = -inf;
    S.defaultBandNames = {};
    S.defaultBandTS = [];

    % Required for averaging: store abs scalogram
    S.absWT = [];
end

function gp = finalise_global_power_max(globalPowerMax)
    if ~isfinite(globalPowerMax)
        gp = 0;
    else
        gp = globalPowerMax + 0.8;
        if gp < 0, gp = 0; end
    end
end

%% ----------------------------- PROCESS STREAM ----------------------------
function [specOut, peaksOut, plotMetaAdd, detailCounter, detailIndexRows] = process_one_stream( ...
    x, FB, N, time_day, condChangeIdx, ...
    fileStem, signalID, sourceTag, dirs, NAMES, SAVE_DPI, ...
    lightStateNumVec, photoperiod_h, ...
    SPEC_TEMPLATE, scaloYTicks, scaloHLinesHours, ...
    BAND_CR, BANDS_UR, detailCacheDir, detailCounter, detailIndexRows, ...
    BAND_COLOUR)

    specOut = SPEC_TEMPLATE;
    peaksOut = {};
    plotMetaAdd = {};

    if isempty(x), return; end

    % Normalise copy for qualitative comparability
    xNorm = x(:);
    if std(xNorm,'omitnan') > 0
        xNorm = (xNorm - mean(xNorm,'omitnan')) ./ std(xNorm,'omitnan');
    else
        xNorm = xNorm * 0;
    end

    % CWT absolute
    coi = [];
    try
        [wt, periods, coi] = cwt(x, 'FilterBank', FB);
    catch
        [wt, periods] = cwt(x, 'FilterBank', FB);
        coi = [];
    end
    periods_hours = hours(periods);

    coi_hours = [];
    if ~isempty(coi)
        try, coi_hours = hours(coi(:));
        catch, coi_hours = [];
        end
    end

    absWT = abs(wt);
    pMax = max(absWT(:));
    p99  = prctile(absWT(:), 99);

    climAbs = [0, max(p99, eps)];
    plotMetaAdd(end+1,:) = {fileStem, signalID, sourceTag, 'Abs', climAbs(1), climAbs(2), pMax, p99, 'CLim based on 99th percentile of abs(wt)'}; %#ok<AGROW>

    % Plot abs scalogram
    fig = [];
    try
        fig = figure('Visible','off');
        pcolor(time_day, periods_hours, absWT); shading interp;
        colormap jet; colorbar;
        caxis(climAbs);
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman','YTick', scaloYTicks);
        xlabel('Time (days)'); ylabel('Period (hr)');
        title(sprintf('Scalogram (Abs) - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');

        fn = sprintf(NAMES.FN_WAV_ABS, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.FigWavAbs, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    % CWT normalised
    try
        [wtN, periodsN] = cwt(xNorm, 'FilterBank', FB);
    catch
        [wtN, periodsN] = cwt(xNorm, 'FilterBank', FB);
    end
    absWTN = abs(wtN);
    p99N  = prctile(absWTN(:), 99);
    climN = [0, max(p99N, eps)];
    plotMetaAdd(end+1,:) = {fileStem, signalID, sourceTag, 'NormZ', climN(1), climN(2), max(absWTN(:)), p99N, 'CWT on z-scored signal; CLim=99th percentile'}; %#ok<AGROW>

    fig = [];
    try
        fig = figure('Visible','off');
        pcolor(time_day, periods_hours, absWTN); shading interp;
        colormap jet; colorbar;
        caxis(climN);
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman','YTick', scaloYTicks);
        xlabel('Time (days)'); ylabel('Period (hr)');
        title(sprintf('Scalogram (NormZ) - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');

        fn = sprintf(NAMES.FN_WAV_NORM, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.FigWavNorm, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    % Bandlines plot (use Abs with CLim)
    fig = [];
    try
        fig = figure('Visible','off');
        pcolor(time_day, periods_hours, absWT); shading interp;
        colormap jet; colorbar;
        caxis(climAbs);
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman','YTick', scaloYTicks);
        hold on;

        for k = 1:numel(condChangeIdx)
            r = condChangeIdx(k) + 1;
            if r >= 1 && r <= numel(time_day)
                xlineT = time_day(r);
                plot([xlineT xlineT], [min(periods_hours) max(periods_hours)], 'w:', 'LineWidth', 1.2);
            end
        end

        yMax = max(periods_hours);
        for hh = scaloHLinesHours(:).'
            if hh <= yMax + 1e-6
                plot([min(time_day) max(time_day)], [hh hh], 'w-', 'LineWidth', 0.7);
            end
        end
        hold off;

        xlabel('Time (days)'); ylabel('Period (hr)');
        title(sprintf('Scalogram + Bands - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');
        fn = sprintf(NAMES.FN_WAV_BANDS, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.FigWavBands, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    % log power and spectra
    powerSpec = abs(wt).^2;
    logPow = log10(powerSpec);

    avgPowerSpectrum = mean(logPow, 2);
    localPowerMax = max(avgPowerSpectrum, [], 'omitnan');

    % Condition spectra derived from lightStateNumVec (0/1)
    condLabels = {'Dark','Light'};
    condSpectra = NaN(size(avgPowerSpectrum,1), 2);
    for j = 1:2
        mask = (double(lightStateNumVec(:)) == (j-1));
        if any(mask)
            condSpectra(:,j) = mean(logPow(:, mask), 2);
        end
        localPowerMax = max(localPowerMax, max(condSpectra(:,j), [], 'omitnan'));
    end

    % Peaks from condition spectra (QC)
    if exist('findpeaks','file') == 2
        for j = 1:2
            spec = condSpectra(:,j);
            if all(~isfinite(spec)), continue; end
            [condPeaks, condLocs, condWidths, condProms] = findpeaks(spec, 'MinPeakProminence', 0);
            if isempty(condPeaks), continue; end

            [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
            sortedLocs   = condLocs(sortIdx);
            sortedWidths = condWidths(sortIdx);
            sortedProms  = condProms(sortIdx);

            sortedPeriods = periods_hours(sortedLocs);
            condParsed = parse_condition_from_signalID(signalID);

            for p = 1:numel(sortedPeaks)
                peaksOut(end+1,:) = { ...
                    fileStem, signalID, sourceTag, condParsed, photoperiod_h, condLabels{j}, 'All', ...
                    p, sortedPeriods(p), sortedPeaks(p), sortedProms(p), sortedWidths(p) ...
                    }; %#ok<AGROW>
            end
        end
    end

    % specOut
    specOut.ok = true;
    specOut.signalID = signalID;
    specOut.sourceTag = sourceTag;
    specOut.periods_hours = periods_hours(:);
    specOut.time_day = time_day(:);
    specOut.logPow = logPow;
    specOut.coi_hours = coi_hours(:);
    specOut.avgPowerSpectrum = avgPowerSpectrum(:);
    specOut.condLabels = condLabels;
    specOut.condSpectra = condSpectra;
    specOut.localPowerMax = localPowerMax;

    % store absWT for averaging
    specOut.absWT = absWT;

    % Ridge overlay on Abs scalogram using bandTS
    fig = [];
    try
        if strcmpi(sourceTag,'Raw')
            bands = [reshape(BAND_CR,1,2); BANDS_UR];
            bandNames = {'CR_20_28','UR_1_3','UR_3_6','UR_6_9','UR_9_12','UR_12_18'};
        else
            bands = BANDS_UR;
            bandNames = {'UR_1_3','UR_3_6','UR_6_9','UR_9_12','UR_12_18'};
        end
        bandTS = compute_band_timeseries(periods_hours, logPow, bands, coi_hours, wt);
        specOut.defaultBandNames = bandNames;
        specOut.defaultBandTS = bandTS;
        cRGB = get_band_colours_rgb(bandNames, BAND_COLOUR);

        fig = figure('Visible','off');
        pcolor(time_day, periods_hours, absWT); shading interp;
        colormap jet; colorbar;
        caxis(climAbs);
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman','YTick', scaloYTicks);
        hold on;

        for b = 1:numel(bandNames)
            rp = bandTS.RidgePeriod(b,:);
            % outline + coloured core
            plot(time_day, rp, '-', 'LineWidth', 3.0, 'Color', [0 0 0]);
            if ~isempty(cRGB) && size(cRGB,1) >= b
                plot(time_day, rp, '-', 'LineWidth', 2.0, 'Color', cRGB(b,:));
            else
                plot(time_day, rp, '-', 'LineWidth', 2.0, 'Color', [1 1 1]);
            end
        end

        hold off;
        xlabel('Time (days)'); ylabel('Period (hr)');
        title(sprintf('Scalogram + Ridge - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');
        fn = sprintf(NAMES.FN_WAV_RIDGE, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.FigWavRidge, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    % Detail cache spectra
    detailCounter = detailCounter + 1;
    shortID = sprintf('C%03d', detailCounter);
    specSheet = ['Spec_' shortID];
    detailIndexRows(end+1,:) = {shortID, fileStem, signalID, sourceTag, specSheet}; %#ok<AGROW>

    specStruct = struct();
    specStruct.shortID = shortID;
    specStruct.fileStem = fileStem;
    specStruct.signalID = signalID;
    specStruct.sourceTag = sourceTag;
    specStruct.periods_hours = periods_hours(:);
    specStruct.globalSpec = avgPowerSpectrum(:);
    specStruct.condLabels = condLabels(:);
    specStruct.condSpectra = condSpectra;

    save(fullfile(detailCacheDir, ['Spec_' shortID '.mat']), 'specStruct');

    clear wt wtN powerSpec logPow absWT absWTN;
end

function condParsed = parse_condition_from_signalID(signalID)
    s = char(string(signalID));
    k = find(s=='_', 1, 'first');
    if isempty(k) || k == 1
        condParsed = s;
    else
        condParsed = s(1:k-1);
    end
end

function safe_close(figH)
    try
        if ~isempty(figH) && isvalid(figH), close(figH); end
    catch
    end
end

%% ----------------------------- BAND TS ----------------------------------
function bandTS = compute_band_timeseries(periods_hours, logPow, bands, coi_hours, wt)
    per = double(periods_hours(:));
    nT  = size(logPow, 2);

    if isvector(bands) && numel(bands)==2
        bands = reshape(bands, 1, 2);
    end
    nB = size(bands, 1);

    bandTS = struct();
    bandTS.BandPower   = NaN(nB, nT);
    bandTS.RidgePeriod = NaN(nB, nT);
    bandTS.RidgePower  = NaN(nB, nT);
    bandTS.RidgePhase  = NaN(nB, nT);
    bandTS.ValidFlag   = false(nB, nT);

    hasWT = (nargin >= 5) && ~isempty(wt) && isequal(size(wt), size(logPow));

    useCOI = ~isempty(coi_hours) && numel(coi_hours) >= nT;
    if useCOI
        coiH = double(coi_hours(:));
        coiH = coiH(1:nT);
    else
        coiH = [];
    end

    for b = 1:nB
        Pmin = bands(b,1);
        Pmax = bands(b,2);

        idx = find(per >= Pmin & per <= Pmax);
        if isempty(idx), continue; end

        sub = logPow(idx, :);

        validT = true(1, nT);
        if useCOI
            validT = isfinite(coiH).' & (coiH.' >= Pmax);
        end

        bp = mean(sub, 1, 'omitnan');
        [mx, k] = max(sub, [], 1, 'includenan');

        ok = validT & isfinite(bp) & isfinite(mx) & isfinite(k);

        bp(~ok) = NaN;
        mx(~ok) = NaN;
        k(~ok)  = NaN;

        bandTS.BandPower(b,:) = bp;
        bandTS.RidgePower(b,:) = mx;

        perBand = per(idx);
        rp = NaN(1, nT);
        rp(ok) = perBand(k(ok));
        bandTS.RidgePeriod(b,:) = rp;

        if hasWT && any(ok)
            phaseRow = NaN(1, nT);
            okIdx = find(ok);
            kOK = k(ok);
            for jj = 1:numel(okIdx)
                try
                    phaseRow(okIdx(jj)) = angle(wt(idx(kOK(jj)), okIdx(jj)));
                catch
                    phaseRow(okIdx(jj)) = NaN;
                end
            end
            bandTS.RidgePhase(b,:) = phaseRow;
        end

        bandTS.ValidFlag(b,:) = ok;
    end
end

function [bandPowRows, ridgePRows, ridgePowRows] = append_band_timeseries_rows_with_valid( ...
    bandPowRows, ridgePRows, ridgePowRows, ...
    fileStem, signalID, sourceTag, ...
    bandNames, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTS)

    nT = numel(time_day);
    nB = numel(bandNames);

    phaseList = {'All'};
    if isfield(phaseMasks,'Light') && isfield(phaseMasks,'Dark')
        phaseList = {'All','Light','Dark'};
    end

    for p = 1:numel(phaseList)
        phaseTag = phaseList{p};
        mask = resize_to_N_logical(phaseMasks.(phaseTag), nT);

        for b = 1:nB
            bn = bandNames{b};

            vBP  = bandTS.BandPower(b,:).';
            vRP  = bandTS.RidgePeriod(b,:).';
            vRpw = bandTS.RidgePower(b,:).';
            vOK  = bandTS.ValidFlag(b,:).';

            for t = 1:nT
                if ~mask(t), continue; end
                bandPowRows(end+1,:)  = {fileStem, signalID, sourceTag, bn, time_day(t), ZT_hr(t), lightStateStrVec{t}, phaseTag, vBP(t), double(vOK(t))}; %#ok<AGROW>
                ridgePRows(end+1,:)   = {fileStem, signalID, sourceTag, bn, time_day(t), ZT_hr(t), lightStateStrVec{t}, phaseTag, vRP(t), double(vOK(t))}; %#ok<AGROW>
                ridgePowRows(end+1,:) = {fileStem, signalID, sourceTag, bn, time_day(t), ZT_hr(t), lightStateStrVec{t}, phaseTag, vRpw(t), double(vOK(t))}; %#ok<AGROW>
            end
        end
    end
end

function m2 = resize_to_N_logical(m, N)
    m = logical(m(:));
    if numel(m) < N
        if isempty(m)
            m2 = false(N,1);
        else
            m2 = [m; repmat(m(end), N-numel(m), 1)];
        end
    elseif numel(m) > N
        m2 = m(1:N);
    else
        m2 = m;
    end
end

function rows = append_period_candidate_rows(fileStem, signalID, sourceTag, hsubModeTag, photoperiod_h, phaseMasks, bandNames, bandTS, minCoverage, minCOIFrac)
    rows = {};
    if isempty(bandNames) || isempty(bandTS), return; end

    nT = size(bandTS.BandPower, 2);
    phaseList = {'All'};
    if isfield(phaseMasks,'Light') && isfield(phaseMasks,'Dark')
        phaseList = {'All','Light','Dark'};
    end

    condParsed = parse_condition_from_signalID(signalID);

    for p = 1:numel(phaseList)
        phaseTag = phaseList{p};
        mask = resize_to_N_logical(phaseMasks.(phaseTag), nT);
        totalN = sum(mask);

        for b = 1:numel(bandNames)
            bn = bandNames{b};
            vOK = logical(bandTS.ValidFlag(b,:).') & mask(:);
            rp  = bandTS.RidgePeriod(b,:).';
            bp  = bandTS.BandPower(b,:).';
            rpw = bandTS.RidgePower(b,:).';

            validN = sum(vOK & isfinite(rp));
            if totalN > 0
                coverage = validN / totalN;
            else
                coverage = NaN;
            end
            coiFrac = coverage; % ValidFlag is the current COI-valid time mask after band/QC filtering.

            medRP = median(rp(vOK), 'omitnan');
            iqrRP = local_iqr_omitnan(rp(vOK));
            meanBP = mean(bp(vOK), 'omitnan');
            sdBP = std(bp(vOK), 0, 'omitnan');
            meanRpw = mean(rpw(vOK), 'omitnan');
            sdRpw = std(rpw(vOK), 0, 'omitnan');

            passQC = isfinite(medRP) && coverage >= minCoverage && coiFrac >= minCOIFrac;
            if passQC
                qcReason = 'Pass';
            elseif ~isfinite(medRP)
                qcReason = 'NoFiniteRidgePeriod';
            elseif coverage < minCoverage
                qcReason = sprintf('LowRidgeCoverage_<%.3g', minCoverage);
            elseif coiFrac < minCOIFrac
                qcReason = sprintf('LowCOIValidFrac_<%.3g', minCOIFrac);
            else
                qcReason = 'FailQC';
            end

            candID = make_candidate_id(fileStem, signalID, sourceTag, hsubModeTag, phaseTag, bn, 1);
            rows(end+1,:) = {fileStem, signalID, condParsed, sourceTag, hsubModeTag, photoperiod_h, phaseTag, bn, ...
                             candID, 1, medRP, iqrRP, meanBP, sdBP, meanRpw, sdRpw, ...
                             coverage, coiFrac, validN, totalN, double(passQC), qcReason}; %#ok<AGROW>
        end
    end
end

function rows = append_ridge_trajectory_rows(fileStem, signalID, sourceTag, hsubModeTag, photoperiod_h, bandNames, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTS)
    rows = {};
    if isempty(bandNames) || isempty(bandTS), return; end

    nT = numel(time_day);
    phaseList = {'All'};
    if isfield(phaseMasks,'Light') && isfield(phaseMasks,'Dark')
        phaseList = {'All','Light','Dark'};
    end
    condParsed = parse_condition_from_signalID(signalID);

    for p = 1:numel(phaseList)
        phaseTag = phaseList{p};
        mask = resize_to_N_logical(phaseMasks.(phaseTag), nT);
        for b = 1:numel(bandNames)
            bn = bandNames{b};
            candID = make_candidate_id(fileStem, signalID, sourceTag, hsubModeTag, phaseTag, bn, 1);
            bp  = bandTS.BandPower(b,:).';
            rp  = bandTS.RidgePeriod(b,:).';
            rpw = bandTS.RidgePower(b,:).';
            ok  = bandTS.ValidFlag(b,:).';
            for t = 1:nT
                if ~mask(t), continue; end
                rows(end+1,:) = {fileStem, signalID, condParsed, sourceTag, hsubModeTag, photoperiod_h, bn, candID, ...
                                 time_day(t), ZT_hr(t), lightStateStrVec{t}, phaseTag, bp(t), rp(t), rpw(t), double(ok(t))}; %#ok<AGROW>
            end
        end
    end
end

function rows = append_ridge_phase_rows(fileStem, signalID, sourceTag, hsubModeTag, photoperiod_h, bandNames, time_day, ZT_hr, lightStateStrVec, phaseMasks, bandTS)
    rows = {};
    if isempty(bandNames) || isempty(bandTS) || ~isfield(bandTS,'RidgePhase'), return; end

    nT = numel(time_day);
    phaseList = {'All'};
    if isfield(phaseMasks,'Light') && isfield(phaseMasks,'Dark')
        phaseList = {'All','Light','Dark'};
    end
    condParsed = parse_condition_from_signalID(signalID);

    for p = 1:numel(phaseList)
        phaseTag = phaseList{p};
        mask = resize_to_N_logical(phaseMasks.(phaseTag), nT);
        for b = 1:numel(bandNames)
            bn = bandNames{b};
            candID = make_candidate_id(fileStem, signalID, sourceTag, hsubModeTag, phaseTag, bn, 1);
            rp  = bandTS.RidgePeriod(b,:).';
            rpw = bandTS.RidgePower(b,:).';
            ph  = bandTS.RidgePhase(b,:).';
            ok  = bandTS.ValidFlag(b,:).';
            for t = 1:nT
                if ~mask(t), continue; end
                rows(end+1,:) = {fileStem, signalID, condParsed, sourceTag, hsubModeTag, photoperiod_h, bn, candID, ...
                                 time_day(t), ZT_hr(t), lightStateStrVec{t}, phaseTag, rp(t), rpw(t), ph(t), double(ok(t))}; %#ok<AGROW>
            end
        end
    end
end

function id = make_candidate_id(fileStem, signalID, sourceTag, hsubModeTag, phaseTag, bandName, rank)
    raw = sprintf('%s__%s__%s__%s__%s__%s__R%d', char(string(fileStem)), char(string(signalID)), ...
                  char(string(sourceTag)), char(string(hsubModeTag)), char(string(phaseTag)), char(string(bandName)), rank);
    id = sanitise_filename(raw);
end

function val = local_iqr_omitnan(x)
    x = x(:);
    x = x(isfinite(x));
    if isempty(x)
        val = NaN;
    else
        try
            val = iqr(x);
        catch
            val = prctile(x,75) - prctile(x,25);
        end
    end
end

function rows = append_band_condition_summary_linearfrac(fileStem, signalID, sourceTag, photoperiod_h, phaseMasks, bandNames, bandTS)
    rows = {};
    nB = numel(bandNames);
    nT = size(bandTS.BandPower, 2);

    phaseList = {'All'};
    if isfield(phaseMasks,'Light') && isfield(phaseMasks,'Dark')
        phaseList = {'All','Light','Dark'};
    end

    for p = 1:numel(phaseList)
        phaseTag = phaseList{p};
        mask = resize_to_N_logical(phaseMasks.(phaseTag), nT);
        if ~any(mask), continue; end

        meanBP_log = NaN(nB,1);
        sdBP_log   = NaN(nB,1);
        meanBP_lin = NaN(nB,1);
        sdBP_lin   = NaN(nB,1);

        meanRP     = NaN(nB,1);
        sdRP       = NaN(nB,1);

        meanRpw_log= NaN(nB,1);
        sdRpw_log  = NaN(nB,1);

        for b = 1:nB
            vBP_log  = bandTS.BandPower(b, mask);
            vRP      = bandTS.RidgePeriod(b, mask);
            vRpw_log = bandTS.RidgePower(b, mask);

            meanBP_log(b) = mean(vBP_log, 'omitnan');
            sdBP_log(b)   = std(vBP_log,  'omitnan');

            vBP_lin = 10 .^ vBP_log;
            meanBP_lin(b) = mean(vBP_lin, 'omitnan');
            sdBP_lin(b)   = std(vBP_lin,  'omitnan');

            meanRP(b) = mean(vRP, 'omitnan');
            sdRP(b)   = std(vRP,  'omitnan');

            meanRpw_log(b) = mean(vRpw_log, 'omitnan');
            sdRpw_log(b)   = std(vRpw_log,  'omitnan');
        end

        denom = sum(meanBP_lin, 'omitnan');
        if ~isfinite(denom) || denom <= 0
            fracLin = NaN(nB,1);
        else
            fracLin = meanBP_lin ./ denom;
        end

        lightStateValStr = 'DerivedFromLightDuration';

        for b = 1:nB
            rows(end+1,:) = { ...
                fileStem, signalID, sourceTag, photoperiod_h, lightStateValStr, phaseTag, bandNames{b}, ...
                meanBP_log(b), sdBP_log(b), ...
                meanBP_lin(b), sdBP_lin(b), ...
                fracLin(b), ...
                meanRP(b), sdRP(b), ...
                meanRpw_log(b), sdRpw_log(b) ...
                }; %#ok<AGROW>
        end
    end
end

%% ----------------------------- PLOTS ------------------------------------
function plot_bandridge_timeseries_figs(fileStem, signalID, sourceTag, time_day, condChangeIdx, bandNames, bandTS, ...
    dirs, NAMES, SAVE_DPI, forceZeroBP_RPOW, forceZeroRP, bandColourRGB)

    legLoc = 'eastoutside';

    % BandPower (log10)
    fig = [];
    try
        fig = figure('Visible','off'); set(fig,'Position',[100 100 1100 420]);
        hold on;
        for b = 1:numel(bandNames)
            if ~isempty(bandColourRGB) && size(bandColourRGB,1) >= b
                plot(time_day, bandTS.BandPower(b,:), 'LineWidth', 1.2, 'Color', bandColourRGB(b,:));
            else
                plot(time_day, bandTS.BandPower(b,:), 'LineWidth', 1.2);
            end
        end
        draw_boundaries(time_day, condChangeIdx);
        hold off;
        xlabel('Time (days)'); ylabel('BandPower (log_{10})');
        title(sprintf('BandPower(t) - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman');
        if forceZeroBP_RPOW
            yl = ylim; ylim([0, max(yl(2), 0.001)]);
        end
        lg = legend(bandNames,'Location',legLoc,'Interpreter','none'); set(lg,'Box','off');
        fn = sprintf(NAMES.FN_BP_COMB, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.BR_BP_Comb, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    % RidgePeriod
    fig = [];
    try
        fig = figure('Visible','off'); set(fig,'Position',[100 100 1100 420]);
        hold on;
        for b = 1:numel(bandNames)
            if ~isempty(bandColourRGB) && size(bandColourRGB,1) >= b
                plot(time_day, bandTS.RidgePeriod(b,:), 'LineWidth', 1.2, 'Color', bandColourRGB(b,:));
            else
                plot(time_day, bandTS.RidgePeriod(b,:), 'LineWidth', 1.2);
            end
        end
        draw_boundaries(time_day, condChangeIdx);
        hold off;
        xlabel('Time (days)'); ylabel('RidgePeriod (hr)');
        title(sprintf('RidgePeriod(t) - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman');
        if forceZeroRP
            yl = ylim; ylim([0, max(yl(2), 0.001)]);
        end
        lg = legend(bandNames,'Location',legLoc,'Interpreter','none'); set(lg,'Box','off');
        fn = sprintf(NAMES.FN_RP_COMB, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.BR_RP_Comb, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    % RidgePower (log10)
    fig = [];
    try
        fig = figure('Visible','off'); set(fig,'Position',[100 100 1100 420]);
        hold on;
        for b = 1:numel(bandNames)
            if ~isempty(bandColourRGB) && size(bandColourRGB,1) >= b
                plot(time_day, bandTS.RidgePower(b,:), 'LineWidth', 1.2, 'Color', bandColourRGB(b,:));
            else
                plot(time_day, bandTS.RidgePower(b,:), 'LineWidth', 1.2);
            end
        end
        draw_boundaries(time_day, condChangeIdx);
        hold off;
        xlabel('Time (days)'); ylabel('RidgePower (log_{10})');
        title(sprintf('RidgePower(t) - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman');
        if forceZeroBP_RPOW
            yl = ylim; ylim([0, max(yl(2), 0.001)]);
        end
        lg = legend(bandNames,'Location',legLoc,'Interpreter','none'); set(lg,'Box','off');
        fn = sprintf(NAMES.FN_RPOW_COMB, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(dirs.BR_RPOW_Comb, fn), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);
end

function draw_boundaries(time_day, condChangeIdx)
    if isempty(condChangeIdx), return; end
    yl = ylim;
    for k = 1:numel(condChangeIdx)
        r = condChangeIdx(k) + 1;
        if r >= 1 && r <= numel(time_day)
            xlineT = time_day(r);
            plot([xlineT xlineT], yl, 'k:', 'LineWidth', 1.0);
        end
    end
end

function plot_power_plots_for_stream(S, fileStem, sourceTag, dirs, NAMES, SAVE_DPI, globalPowerMax, globalYMaxHr)
    signalID = S.signalID;
    periods_hours = S.periods_hours;
    avgSpec = S.avgPowerSpectrum;
    condSpectra = S.condSpectra;
    condLabels = S.condLabels;

    if strcmpi(sourceTag,'Raw')
        outGlobDir = dirs.FigPwr_Raw_Glob;
        outCondDir = dirs.FigPwr_Raw_Cond;
    else
        outGlobDir = dirs.FigPwr_H_Glob;
        outCondDir = dirs.FigPwr_H_Cond;
    end

    fig = [];
    try
        fig = figure('Visible','off');
        plot(avgSpec, periods_hours, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})'); ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - %s | %s | %s', signalID, fileStem, sourceTag), 'Interpreter','none');
        set(gca, 'TickDir','out', 'FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]); ylim([0, globalYMaxHr]);
        fnG = sprintf(NAMES.FN_PWR_GLOB, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag);
        print(fig, fullfile(outGlobDir, fnG), '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);

    for j = 1:numel(condLabels)
        spec = condSpectra(:, j);
        if all(~isfinite(spec)), continue; end
        fig = [];
        try
            fig = figure('Visible','off');
            plot(spec, periods_hours, '-k', 'LineWidth', 1.5);
            xlabel('Power (log_{10})'); ylabel('Period (hr)');
            title(sprintf('Power Spectrum - %s | %s | %s | %s', signalID, fileStem, sourceTag, condLabels{j}), 'Interpreter','none');
            set(gca, 'TickDir','out', 'FontName','Times New Roman','box','off');
            xlim([0, globalPowerMax]); ylim([0, globalYMaxHr]);
            condTag = sanitise_filename(condLabels{j});
            fnC = sprintf(NAMES.FN_PWR_COND, sanitise_filename(fileStem), sanitise_filename(signalID), sourceTag, condTag);
            print(fig, fullfile(outCondDir, fnC), '-djpeg', sprintf('-r%d', SAVE_DPI));
        catch
        end
        safe_close(fig);
    end
end

function row = compute_band_summary_row(fileStem, signalID, sourceTag, photoperiod_h, periods_hours, avgSpec, BAND_CR, BANDS_UR)
    row = {};
    if isempty(periods_hours) || isempty(avgSpec), return; end
    per = double(periods_hours(:));
    sp  = double(avgSpec(:));

    meanCR = NaN;
    urVals = nan(1, size(BANDS_UR,1));

    if strcmpi(sourceTag, 'Raw')
        mask = per >= BAND_CR(1) & per <= min(BAND_CR(2), max(per));
        if any(mask), meanCR = mean(sp(mask), 'omitnan'); end
    elseif strcmpi(sourceTag, 'HSub')
        for b = 1:size(BANDS_UR,1)
            lo = BANDS_UR(b,1); hi = BANDS_UR(b,2);
            mask = per >= lo & per <= hi;
            if any(mask), urVals(b) = mean(sp(mask), 'omitnan'); end
        end
    end

    row = {fileStem, signalID, sourceTag, photoperiod_h, meanCR, urVals(1), urVals(2), urVals(3), urVals(4), urVals(5)};
end

%% ----------------------------- WORKBOOKS --------------------------------
function write_summary_workbook(summaryXLSX, InputsTable, PeaksTable, BandSummaryTable, BandConditionSummaryTable, PeriodCandidatesTable, periodMinMin, periodMaxMin, topN, runType, doPhase, primaryHSubResidual, minRidgeCoverage, minCOIValidFrac)
    readme = {
        'Wavelet_analysis summary'
        ''
        sprintf('RunType: %s', char(string(runType)))
        sprintf('PhaseSplitEnabled: %s', string(doPhase))
        ''
        'Core wavelet method'
        sprintf('- PeriodLimits = [%d, %d] minutes (%.1f to %.1f hours)', periodMinMin, periodMaxMin, periodMinMin/60, periodMaxMin/60)
        '- Light/Dark derived ONLY from Light duration (h) with ZT0=lights-on'
        '- BandConditionSummary: includes log10 and linear band-power summaries; fractions are linear'
        '- v12: Raw exports CR and UR bands; HSub exports UR bands for validation only'
        sprintf('- Primary HSub validation residual when v12 manifest is available: %s', char(string(primaryHSubResidual)))
        sprintf('- Period candidate PassQC requires RidgeCoverageFrac >= %.3g and COIValidFrac >= %.3g', minRidgeCoverage, minCOIValidFrac)
        };
    writecell(readme(:), summaryXLSX, 'Sheet', 'README');
    writetable(InputsTable, summaryXLSX, 'Sheet', 'Inputs');

    if ~isempty(PeaksTable)
        PeaksTable2 = sortrows(PeaksTable, 'PeakValue_log10', 'descend');
        writetable(PeaksTable2, summaryXLSX, 'Sheet', 'Peaks');
    else
        writecell({'No peaks extracted.'}, summaryXLSX, 'Sheet', 'Peaks');
    end

    if ~isempty(BandSummaryTable)
        writetable(BandSummaryTable, summaryXLSX, 'Sheet', 'BandSummary');
    else
        writecell({'No BandSummary rows produced.'}, summaryXLSX, 'Sheet', 'BandSummary');
    end

    if ~isempty(BandConditionSummaryTable)
        writetable(BandConditionSummaryTable, summaryXLSX, 'Sheet', 'BandConditionSummary');
    else
        writecell({'No BandConditionSummary rows produced.'}, summaryXLSX, 'Sheet', 'BandConditionSummary');
    end

    if ~isempty(PeriodCandidatesTable)
        writetable(PeriodCandidatesTable, summaryXLSX, 'Sheet', 'PeriodCandidates');
    else
        writecell({'No PeriodCandidates rows produced.'}, summaryXLSX, 'Sheet', 'PeriodCandidates');
    end
end

function write_detail_workbook(detailXLSX, DetailIndexTable, detailCacheDir)
    writecell({'Detail spectra workbook. Uses MAT cache files in DetailCache.'}, detailXLSX, 'Sheet', 'README');

    if ~isempty(DetailIndexTable)
        writetable(DetailIndexTable, detailXLSX, 'Sheet', 'Index');
    else
        writecell({'No detail spectra stored.'}, detailXLSX, 'Sheet', 'Index');
        return;
    end

    for i = 1:height(DetailIndexTable)
        shortID = char(string(DetailIndexTable.ShortID(i)));
        sheet   = char(string(DetailIndexTable.SpectrumSheet(i)));
        matPath = fullfile(detailCacheDir, ['Spec_' shortID '.mat']);
        if ~isfile(matPath), continue; end

        S = load(matPath, 'specStruct');
        specStruct = S.specStruct;

        per = specStruct.periods_hours(:);
        g   = specStruct.globalSpec(:);
        TT = table(per, g, 'VariableNames', {'Period_hr','GlobalPower_log10'});

        condLabels = specStruct.condLabels;
        condSpectra = specStruct.condSpectra;

        if ~isempty(condLabels) && ~isempty(condSpectra)
            for j = 1:numel(condLabels)
                nm = matlab.lang.makeValidName(['Cond_' char(string(condLabels{j})) '_Power_log10']);
                TT.(nm) = condSpectra(:, j);
            end
        end

        if strlength(string(sheet)) > 31
            sheet = char(extractBefore(string(sheet), 32));
        end
        writetable(TT, detailXLSX, 'Sheet', sheet);
    end
end

function write_band_timeseries_workbook(outXLSX, bandPowRows, ridgePRows, ridgePowRows, hdr, RidgeTrajectoryTable, RidgePhaseTable)
    writecell({'WP_BandTimeSeries: long-format band and ridge time-series. Includes ZT and ValidFlag.'}, outXLSX, 'Sheet', 'README');

    if isempty(bandPowRows)
        writecell({'No BandPower rows.'}, outXLSX, 'Sheet', 'BandPower_Long');
    else
        T = cell2table(bandPowRows, 'VariableNames', hdr);
        writetable(T, outXLSX, 'Sheet', 'BandPower_Long');
    end

    if isempty(ridgePRows)
        writecell({'No RidgePeriod rows.'}, outXLSX, 'Sheet', 'RidgePeriod_Long');
    else
        T = cell2table(ridgePRows, 'VariableNames', hdr);
        writetable(T, outXLSX, 'Sheet', 'RidgePeriod_Long');
    end

    if isempty(ridgePowRows)
        writecell({'No RidgePower rows.'}, outXLSX, 'Sheet', 'RidgePower_Long');
    else
        T = cell2table(ridgePowRows, 'VariableNames', hdr);
        writetable(T, outXLSX, 'Sheet', 'RidgePower_Long');
    end

    if ~isempty(RidgeTrajectoryTable)
        writetable(RidgeTrajectoryTable, outXLSX, 'Sheet', 'RidgeTrajectory_Long');
    else
        writecell({'No RidgeTrajectory rows.'}, outXLSX, 'Sheet', 'RidgeTrajectory_Long');
    end

    if ~isempty(RidgePhaseTable)
        writetable(RidgePhaseTable, outXLSX, 'Sheet', 'RidgePhase_Long');
    else
        writecell({'No RidgePhase rows.'}, outXLSX, 'Sheet', 'RidgePhase_Long');
    end
end

%% ========================================================================
% Averaged Abs scalograms (per file): compute only (no MAT), save JPEGs with
% Raw/HSub coupled CLim per group.
% ========================================================================

function avgOut = compute_avg_abs_scalograms_perfile(streamStore, condMap, fileStem, sourceTag)
    avgOut = struct();
    avgOut.ok = false;
    avgOut.fileStem = char(string(fileStem));
    avgOut.sourceTag = char(string(sourceTag));
    avgOut.groups = struct(); % groupLabel -> struct with AvgAbs, Counts, periods, time, nSignals
    avgOut.periods_hours = [];
    avgOut.time_day = [];

    % Eligible signals
    idxOk = [];
    for i = 1:numel(streamStore)
        S = streamStore(i);
        if ~isfield(S,'ok') || ~S.ok, continue; end
        if isempty(S.absWT) || isempty(S.periods_hours) || isempty(S.time_day), continue; end
        idxOk(end+1) = i; %#ok<AGROW>
    end
    if isempty(idxOk)
        return;
    end

    % Enforce common grid within file
    ref = streamStore(idxOk(1));
    refP = ref.periods_hours(:);
    refT = ref.time_day(:);
    nP = numel(refP);
    nT = numel(refT);

    keep = false(size(idxOk));
    for k = 1:numel(idxOk)
        S = streamStore(idxOk(k));
        keep(k) = (numel(S.periods_hours)==nP) && (numel(S.time_day)==nT) && all(abs(double(S.periods_hours(:))-double(refP))<1e-12);
    end
    idxOk = idxOk(keep);
    if isempty(idxOk)
        return;
    end

    avgOut.periods_hours = refP;
    avgOut.time_day = refT;

    % Group labels from condMap
    labs = {};
    if ~isempty(condMap) && istable(condMap)
        labs = unique(cellstr(condMap.ConditionLabel));
    end
    if isempty(labs)
        labs = {'Unassigned'};
    end
    if ~any(strcmp(labs,'Unassigned'))
        labs = [labs; {'Unassigned'}];
    end

    groupList = [{'Overall'}; labs(:)];

    for g = 1:numel(groupList)
        glab = groupList{g};

        sigIdx = [];
        for k = 1:numel(idxOk)
            S = streamStore(idxOk(k));
            lab = get_condition_label(condMap, S.signalID);
            if strcmpi(glab,'Overall')
                sigIdx(end+1) = idxOk(k); %#ok<AGROW>
            else
                if strcmp(lab, glab)
                    sigIdx(end+1) = idxOk(k); %#ok<AGROW>
                end
            end
        end

        if isempty(sigIdx)
            continue;
        end

        % Average Abs scalograms (finite only)
        Sum = zeros(nP,nT);
        Cnt = zeros(nP,nT);

        for k = 1:numel(sigIdx)
            S = streamStore(sigIdx(k));
            A = S.absWT;
            m = isfinite(A);
            A(~m) = 0;
            Sum = Sum + A;
            Cnt = Cnt + double(m);
        end

        Avg = Sum ./ max(Cnt,1);
        Avg(Cnt==0) = NaN;

        out = struct();
        out.fileStem = char(string(fileStem));
        out.sourceTag = char(string(sourceTag));
        out.groupLabel = glab;
        out.periods_hours = refP;
        out.time_day = refT;
        out.AvgAbs = Avg;
        out.Counts = Cnt;
        out.nSignals = numel(sigIdx);

        avgOut.groups.(matlab.lang.makeValidName(glab)) = out;
    end

    avgOut.ok = ~isempty(fieldnames(avgOut.groups));
end

function save_perfile_avgabs_jpegs(avgRaw, avgH, dirs, NAMES, SAVE_DPI)
    gAll = union(fieldnames_safe(avgRaw,'groups'), fieldnames_safe(avgH,'groups'));
    if isempty(gAll), return; end

    for i = 1:numel(gAll)
        gName = gAll{i};

        dispLabel = gName;
        if isfield(avgRaw,'groups') && isfield(avgRaw.groups,gName)
            dispLabel = avgRaw.groups.(gName).groupLabel;
        elseif isfield(avgH,'groups') && isfield(avgH.groups,gName)
            dispLabel = avgH.groups.(gName).groupLabel;
        end

        haveR = isfield(avgRaw,'groups') && isfield(avgRaw.groups,gName);
        haveH = isfield(avgH,'groups')   && isfield(avgH.groups,gName);

        climMax = compute_coupled_clim_max( ...
            get_avgabs_if_present(avgRaw, gName), ...
            get_avgabs_if_present(avgH, gName) );

        if haveR
            outDir = dirs.FigWavAvgRaw;
            ensure_dir(outDir);
            fn = sprintf(NAMES.FN_AVG_ABS, sanitise_filename(avgRaw.fileStem), 'Raw', sanitise_filename(dispLabel));
            outPath = fullfile(outDir, fn);
            plot_avg_abs_scalogram_jpeg(avgRaw.groups.(gName).time_day, avgRaw.groups.(gName).periods_hours, avgRaw.groups.(gName).AvgAbs, ...
                climMax, sprintf('AVG Abs scalogram - %s | Raw | %s', avgRaw.fileStem, dispLabel), outPath, SAVE_DPI, []);
        end

        if haveH
            outDir = dirs.FigWavAvgHSub;
            ensure_dir(outDir);
            fn = sprintf(NAMES.FN_AVG_ABS, sanitise_filename(avgH.fileStem), 'HSub', sanitise_filename(dispLabel));
            outPath = fullfile(outDir, fn);
            plot_avg_abs_scalogram_jpeg(avgH.groups.(gName).time_day, avgH.groups.(gName).periods_hours, avgH.groups.(gName).AvgAbs, ...
                climMax, sprintf('AVG Abs scalogram - %s | HSub | %s', avgH.fileStem, dispLabel), outPath, SAVE_DPI, []);
        end
    end
end

function A = get_avgabs_if_present(avgOut, gName)
    A = [];
    try
        if isfield(avgOut,'groups') && isfield(avgOut.groups, gName)
            A = avgOut.groups.(gName).AvgAbs;
        end
    catch
        A = [];
    end
end

function mx = compute_coupled_clim_max(Araw, Ahsub)
    v = [];
    if ~isempty(Araw)
        vr = Araw(:); vr = vr(isfinite(vr));
        v = [v; vr]; %#ok<AGROW>
    end
    if ~isempty(Ahsub)
        vh = Ahsub(:); vh = vh(isfinite(vh));
        v = [v; vh]; %#ok<AGROW>
    end
    if isempty(v)
        mx = 1;
        return;
    end
    mx = prctile(v, 99);
    if ~isfinite(mx) || mx <= 0
        mx = max(v);
        if ~isfinite(mx) || mx <= 0, mx = 1; end
    end
    mx = max(mx, eps);
end

function plot_avg_abs_scalogram_jpeg(time_day, periods_hours, AvgAbs, climMax, ttl, outPath, SAVE_DPI, boundaryTimes)
    fig = [];
    try
        fig = figure('Visible','off');
        pcolor(time_day, periods_hours, AvgAbs); shading interp;
        colormap jet; colorbar;
        caxis([0 climMax]);
        set(gca,'Box','off','TickDir','out','FontName','Times New Roman','YTick', 0:4:26);
        xlabel('Time (days)'); ylabel('Period (hr)');
        title(ttl, 'Interpreter','none');

        if ~isempty(boundaryTimes)
            hold on;
            yMin = min(periods_hours); yMax = max(periods_hours);
            for k = 1:numel(boundaryTimes)
                xB = boundaryTimes(k);
                plot([xB xB], [yMin yMax], 'w:', 'LineWidth', 1.2);
            end
            hold off;
        end

        print(fig, outPath, '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
    safe_close(fig);
end

function lab = get_condition_label(condMap, signalID)
    lab = 'Unassigned';
    if isempty(condMap) || ~istable(condMap), return; end
    try
        m = strcmp(condMap.SignalID, string(signalID));
        if any(m)
            lab = char(condMap.ConditionLabel(find(m,1,'first')));
        end
    catch
        lab = 'Unassigned';
    end
end

function f = fieldnames_safe(S, fieldName)
    f = {};
    try
        if ~isempty(S) && isstruct(S) && isfield(S,fieldName) && isstruct(S.(fieldName))
            f = fieldnames(S.(fieldName));
        end
    catch
        f = {};
    end
end

%% ========================================================================
% Cross-file pooled accumulators (store sums/counts only, write combined JPEGs)
% ========================================================================

function P = init_crosspool()
    P = struct();
    P.Raw = struct();   % photKey -> struct of groups accumulators + reference grid
    P.HSub = struct();
end

function P = crosspool_add_file(P, avgOut, photoperiod_h, sourceTag)
    if isempty(avgOut) || ~isfield(avgOut,'ok') || ~avgOut.ok
        return;
    end

    photKey = photoperiod_key(photoperiod_h); % e.g., L12
    if strcmpi(sourceTag,'Raw')
        P.Raw = crosspool_add_to_bucket(P.Raw, photKey, avgOut);
    else
        P.HSub = crosspool_add_to_bucket(P.HSub, photKey, avgOut);
    end
end

function bucket = crosspool_add_to_bucket(bucket, key, avgOut)
    keyF = matlab.lang.makeValidName(char(string(key)));
    if ~isfield(bucket, keyF)
        bucket.(keyF) = init_crosspool_bucket(avgOut);
    end

    B = bucket.(keyF);

    if ~is_same_grid(B.periods_hours, B.time_day, avgOut.periods_hours, avgOut.time_day)
        bucket.(keyF) = B;
        return;
    end

    gNames = fieldnames(avgOut.groups);
    for i = 1:numel(gNames)
        gn = gNames{i};
        G = avgOut.groups.(gn);

        if ~isfield(B.groups, gn)
            B.groups.(gn) = init_group_accumulator(size(G.AvgAbs));
        end
        A = B.groups.(gn);

        sumAdd = G.AvgAbs .* G.Counts;
        sumAdd(~isfinite(sumAdd)) = 0;

        A.Sum = A.Sum + sumAdd;
        A.Cnt = A.Cnt + G.Counts;
        A.nSignals = A.nSignals + G.nSignals;
        A.nFiles = A.nFiles + 1;

        B.groups.(gn) = A;
    end

    bucket.(keyF) = B;
end

function B = init_crosspool_bucket(avgOut)
    B = struct();
    B.periods_hours = avgOut.periods_hours(:);
    B.time_day = avgOut.time_day(:);
    B.groups = struct();
end

function A = init_group_accumulator(sz)
    A = struct();
    A.Sum = zeros(sz);
    A.Cnt = zeros(sz);
    A.nSignals = 0;
    A.nFiles = 0;
end

function tf = is_same_grid(P1, T1, P2, T2)
    tf = false;
    if isempty(P1) || isempty(T1) || isempty(P2) || isempty(T2), return; end
    if numel(P1) ~= numel(P2), return; end
    if numel(T1) ~= numel(T2), return; end
    if any(abs(double(P1(:)) - double(P2(:))) > 1e-12), return; end
    if any(abs(double(T1(:)) - double(T2(:))) > 1e-12), return; end
    tf = true;
end

function key = photoperiod_key(photoperiod_h)
    if ~isfinite(photoperiod_h)
        key = 'LNaN';
        return;
    end
    L = round(photoperiod_h*10)/10;
    if abs(L - round(L)) < 1e-6
        key = sprintf('L%d', round(L));
    else
        key = sprintf('L%0.1f', L);
        key = strrep(key,'.','p');
    end
end

%% ========================================================================
% Cross-file combined writer (single concatenated scalogram per group)
% Output:
%   parentOut/AcrossFile_Averages/Combined/(Raw|HSub)/  as JPEG only
% ========================================================================

function write_crossfile_combined_abs_jpegs(rootCombinedDir, P, SAVE_DPI)
    if isempty(P) || ~isstruct(P), return; end

    outRaw  = fullfile(rootCombinedDir, 'Raw');  ensure_dir(outRaw);
    outHSub = fullfile(rootCombinedDir, 'HSub'); ensure_dir(outHSub);

    [combRaw, rawMeta]  = build_combined_concat_by_photoperiod(P.Raw);
    [combH,   hsubMeta] = build_combined_concat_by_photoperiod(P.HSub);

    gAll = union(fieldnames_safe(combRaw,'groups'), fieldnames_safe(combH,'groups'));
    if isempty(gAll), return; end

    for i = 1:numel(gAll)
        gn = gAll{i};
        haveR = isfield(combRaw,'groups') && isfield(combRaw.groups, gn);
        haveH = isfield(combH,'groups')   && isfield(combH.groups, gn);
        if ~haveR && ~haveH, continue; end

        dispLabel = gn;
        if haveR, dispLabel = combRaw.groups.(gn).groupLabel; end
        if ~haveR && haveH, dispLabel = combH.groups.(gn).groupLabel; end

        climMax = compute_coupled_clim_max( ...
            get_avgabs_if_present(combRaw, gn), ...
            get_avgabs_if_present(combH, gn) );

        bTimes = [];
        if isfield(rawMeta,'boundaryTimes') && ~isempty(rawMeta.boundaryTimes)
            bTimes = rawMeta.boundaryTimes;
        elseif isfield(hsubMeta,'boundaryTimes') && ~isempty(hsubMeta.boundaryTimes)
            bTimes = hsubMeta.boundaryTimes;
        end

        if haveR
            outPath = fullfile(outRaw, sprintf('COMBINED_AVGABS_%s.jpg', sanitise_filename(dispLabel)));
            plot_avg_abs_scalogram_jpeg(combRaw.time_day, combRaw.periods_hours, combRaw.groups.(gn).AvgAbs, ...
                climMax, sprintf('Combined AVG Abs scalogram | Raw | %s', dispLabel), outPath, SAVE_DPI, bTimes);
        end

        if haveH
            outPath = fullfile(outHSub, sprintf('COMBINED_AVGABS_%s.jpg', sanitise_filename(dispLabel)));
            plot_avg_abs_scalogram_jpeg(combH.time_day, combH.periods_hours, combH.groups.(gn).AvgAbs, ...
                climMax, sprintf('Combined AVG Abs scalogram | HSub | %s', dispLabel), outPath, SAVE_DPI, bTimes);
        end
    end
end

function [comb, meta] = build_combined_concat_by_photoperiod(bucketStruct)
    comb = struct();
    comb.ok = false;
    comb.groups = struct();
    comb.periods_hours = [];
    comb.time_day = [];

    meta = struct();
    meta.photKeysSorted = {};
    meta.boundaryTimes = [];

    if isempty(bucketStruct) || ~isstruct(bucketStruct)
        return;
    end

    keys = fieldnames(bucketStruct);
    if isempty(keys), return; end

    [keysSorted, keysPretty] = sort_photoperiod_keys(keys);
    meta.photKeysSorted = keysPretty;

    refP = [];
    for i = 1:numel(keysSorted)
        B = bucketStruct.(keysSorted{i});
        if isfield(B,'periods_hours') && ~isempty(B.periods_hours)
            refP = B.periods_hours(:);
            break;
        end
    end
    if isempty(refP), return; end
    comb.periods_hours = refP;

    allGroups = {};
    for i = 1:numel(keysSorted)
        B = bucketStruct.(keysSorted{i});
        if isfield(B,'groups') && isstruct(B.groups)
            allGroups = union(allGroups, fieldnames(B.groups));
        end
    end
    if isempty(allGroups), return; end

    timeBlocks = cell(numel(keysSorted),1);
    keepKey = false(numel(keysSorted),1);

    for i = 1:numel(keysSorted)
        B = bucketStruct.(keysSorted{i});
        if ~isfield(B,'periods_hours') || ~isfield(B,'time_day') || isempty(B.periods_hours) || isempty(B.time_day)
            continue;
        end
        if numel(B.periods_hours(:)) ~= numel(refP) || any(abs(double(B.periods_hours(:))-double(refP))>1e-12)
            continue;
        end
        timeBlocks{i} = B.time_day(:);
        keepKey(i) = true;
    end

    keysSorted = keysSorted(keepKey);
    keysPretty = keysPretty(keepKey);
    timeBlocks = timeBlocks(keepKey);

    if isempty(keysSorted), return; end

    dt = [];
    try
        t0 = timeBlocks{1};
        if numel(t0) > 1
            dt = median(diff(t0), 'omitnan');
        end
    catch
        dt = [];
    end
    if isempty(dt) || ~isfinite(dt) || dt <= 0
        dt = 1/24;
    end

    tCat = [];
    boundaryTimes = [];
    offset = 0;

    for i = 1:numel(timeBlocks)
        t = timeBlocks{i};
        if isempty(t), continue; end

        if i == 1
            tAdj = t - t(1);
            offset = tAdj(end) + dt;
        else
            boundaryTimes(end+1) = offset; %#ok<AGROW>
            tAdj = (t - t(1)) + offset;
            offset = tAdj(end) + dt;
        end

        tCat = [tCat; tAdj]; %#ok<AGROW>
    end

    comb.time_day = tCat;
    meta.boundaryTimes = boundaryTimes;
    meta.photKeysSorted = keysPretty;

    for g = 1:numel(allGroups)
        gn = allGroups{g};
        Ablocks = {};

        for i = 1:numel(keysSorted)
            B = bucketStruct.(keysSorted{i});
            if ~isfield(B,'groups') || ~isfield(B.groups, gn), continue; end
            G = B.groups.(gn);
            Avg = G.Sum ./ max(G.Cnt, 1);
            Avg(G.Cnt==0) = NaN;
            Ablocks{end+1} = Avg; %#ok<AGROW>
        end

        if isempty(Ablocks), continue; end

        Acat = [];
        for i = 1:numel(Ablocks)
            Acat = [Acat, Ablocks{i}]; %#ok<AGROW>
        end

        out = struct();
        out.groupLabel = gn;
        out.AvgAbs = Acat;

        comb.groups.(gn) = out;
    end

    comb.ok = ~isempty(fieldnames(comb.groups));
end

function [keysSorted, keysPretty] = sort_photoperiod_keys(keys)
    keysPretty = keys(:);
    Lval = NaN(numel(keys),1);
    for i = 1:numel(keys)
        s = char(keys{i});
        s2 = s;
        if startsWith(s2,'L')
            s2 = s2(2:end);
        end
        s2 = strrep(s2, 'p', '.');
        v = str2double(s2);
        if isfinite(v)
            Lval(i) = v;
        end
    end

    have = isfinite(Lval);
    if any(have)
        [~, ord] = sortrows([Lval(:), double(~have(:))], [1 2]);
        keysSorted = keys(ord);
        keysPretty = keysPretty(ord);
    else
        [keysPretty, ord] = sort(keysPretty);
        keysSorted = keys(ord);
    end
end

%% ========================================================================
% Paul Tol helper functions + sidecar map writer
% ========================================================================

function rgb = hex2rgb_array(hexList)
    n = numel(hexList);
    rgb = NaN(n,3);
    for i = 1:n
        rgb(i,:) = hex2rgb_one(hexList{i});
    end
end

function rgb = hex2rgb_one(hx)
    s = char(string(hx));
    s = strrep(s,'#','');
    if numel(s) ~= 6
        rgb = [0 0 0];
        return;
    end
    r = hex2dec(s(1:2));
    g = hex2dec(s(3:4));
    b = hex2dec(s(5:6));
    rgb = [r g b] / 255;
end

function C = get_band_colours_rgb(bandNames, BAND_COLOUR)
    n = numel(bandNames);
    C = NaN(n,3);
    for i = 1:n
        bn = char(string(bandNames{i}));
        fld = matlab.lang.makeValidName(bn);
        if isfield(BAND_COLOUR, fld)
            C(i,:) = BAND_COLOUR.(fld);
        else
            C(i,:) = [0 0 0];
        end
    end
end

function write_band_colour_map_xlsx(outXLSX, BAND_COLOUR, TOL)
    if exist(outXLSX,'file'), delete(outXLSX); end

    writecell({ ...
        'Band colour map (Paul Tol Bright palette)'; ...
        'Intended as a shared input for downstream scripts (e.g., Script 2).'; ...
        'Colours are stored as Hex and RGB (0-1).'; ...
        }, outXLSX, 'Sheet', 'README');

    pal = cell(numel(TOL.BrightHex), 5);
    for i = 1:numel(TOL.BrightHex)
        pal{i,1} = i;
        pal{i,2} = TOL.BrightHex{i};
        pal{i,3} = TOL.BrightRGB(i,1);
        pal{i,4} = TOL.BrightRGB(i,2);
        pal{i,5} = TOL.BrightRGB(i,3);
    end
    Tpal = cell2table(pal, 'VariableNames', {'Index','Hex','R','G','B'});
    writetable(Tpal, outXLSX, 'Sheet', 'Tol_Bright');

    bNames = fieldnames(BAND_COLOUR);
    rows = cell(numel(bNames), 5);
    for i = 1:numel(bNames)
        bn = bNames{i};
        rgb01 = BAND_COLOUR.(bn);
        rows{i,1} = bn;
        rows{i,2} = rgb01(1);
        rows{i,3} = rgb01(2);
        rows{i,4} = rgb01(3);
        rows{i,5} = rgb2hex(rgb01);
    end
    Tband = cell2table(rows, 'VariableNames', {'BandName','R','G','B','Hex'});
    writetable(Tband, outXLSX, 'Sheet', 'BandColours');
end

function hx = rgb2hex(rgb01)
    rgb01 = max(0,min(1,double(rgb01)));
    v = round(rgb01 * 255);
    hx = sprintf('#%02X%02X%02X', v(1), v(2), v(3));
end