%% Plot_ValidatedUR_PhaseEventHistograms_v1_2.m
% -------------------------------------------------------------------------
% Purpose
%   Publication-style visualisation of validated ultradian ridge-phase events.
%
%   This script uses the validated-Raw pipeline outputs to generate figures
%   inspired by onset-distribution plots, but adapted to the present
%   wavelet-ridge method.
%
%   It does NOT detect behavioural onsets. Instead, it identifies repeated
%   phase-marker events from validated raw ultradian ridge-following phase:
%
%       Event = upward 0-radian phase crossing along a validated raw ridge
%
%   These events are then plotted across ZT and around LD/DL or projected
%   LL anchors.
%
% Inputs
%   1) AcrossPhotoperiod_Input folder containing:
%        WP_TS__*.mat
%      Each file should contain:
%        pkgTS.tables.RidgePhase_Long
%
%   2) HSubSupported_PeriodMap.mat from script C:
%        RawVsSelectiveHSub_PeriodValidation/HSubSupported_PeriodMap.mat
%
% Outputs
%   Plot_ValidatedUR_PhaseEvents/
%       ValidatedUR_PhaseEvents_Output.xlsx
%       Figures/
%           ZT_Histograms/
%           ZT_Rasters/
%           TransitionAligned/
%           Combined/
%
% Notes
%   - Raw signal only.
%   - Uses CarryForward validated candidates only.
%   - LL is plotted using projected previous-L22 anchors only for reference.
%   - These figures should be described as "validated ridge-phase events",
%     not behavioural activity onsets.
%
% Developed for Isaiah J. Ting photoperiod mouse-behaviour pipeline
% -------------------------------------------------------------------------

clear; clc; close all;

%% --------------------------- USER SETTINGS ------------------------------

SELECTED_BANDS = ["UR_1_3","UR_3_6"];  % recommended starting point
SELECTED_PHOTOPERIODS = [12 22 24];    % 24 = LL

% Histogram settings
ZT_BIN_WIDTH_H = 0.25;                 % 15 min bins
REL_BIN_WIDTH_H = 0.25;                % transition-aligned bins
REL_WINDOW_H_DEFAULT = 6;              % +/- window around anchors

% LL projected phase reference
LL_PROJECTED_REFERENCE_PHOTOPERIOD_H = 22;

% Validity filters
REQUIRE_COI_VALID = true;
REQUIRE_VALID_FLAG = true;
MIN_EVENTS_PER_GROUP_FOR_RASTER = 1;

% Output style
FIG_DPI = 600;
FONT_NAME = 'Times New Roman';
FONT_SIZE_AX = 18;
FONT_SIZE_LABEL = 22;
FONT_SIZE_TITLE = 24;
LINE_WIDTH_AX = 1.5;

% Plot options
SAVE_TIFF = true;
SAVE_JPEG = true;
MAKE_RASTERS = true;
MAKE_TRANSITION_ALIGNED = true;
MAKE_COMBINED_ZT_FIGURES = true;

%% ------------------------- SELECT INPUT FOLDER --------------------------

handoffDir = uigetdir(pwd, 'Select AcrossPhotoperiod_Input folder containing WP_TS__*.mat');
if isequal(handoffDir, 0)
    error('No handoff folder selected. Script stopped.');
end

defaultValMap = fullfile(handoffDir, 'RawVsSelectiveHSub_PeriodValidation', ...
    'HSubSupported_PeriodMap.mat');

if exist(defaultValMap, 'file')
    validationMapFile = defaultValMap;
else
    [vf, vp] = uigetfile({'*.mat','MAT-files (*.mat)'}, ...
        'Select HSubSupported_PeriodMap.mat');
    if isequal(vf,0)
        error('No validation map selected. Script stopped.');
    end
    validationMapFile = fullfile(vp, vf);
end

outDir = fullfile(handoffDir, 'Plot_ValidatedUR_PhaseEvents');
figDir = fullfile(outDir, 'Figures');
ztHistDir = fullfile(figDir, 'ZT_Histograms');
rasterDir = fullfile(figDir, 'ZT_Rasters');
relDir = fullfile(figDir, 'TransitionAligned');
combDir = fullfile(figDir, 'Combined');

make_dir(outDir);
make_dir(figDir);
make_dir(ztHistDir);
make_dir(rasterDir);
make_dir(relDir);
make_dir(combDir);

logFile = fullfile(outDir, 'Plot_ValidatedUR_PhaseEvents_Log.txt');
diary(logFile);
fprintf('Plot_ValidatedUR_PhaseEventHistograms_v1_2 started at %s\n', datestr(now));
fprintf('Handoff folder: %s\n', handoffDir);
fprintf('Validation map: %s\n\n', validationMapFile);

%% ------------------------------ STYLE -----------------------------------

set(groot, 'defaultAxesFontName', FONT_NAME);
set(groot, 'defaultTextFontName', FONT_NAME);
set(groot, 'defaultAxesFontSize', FONT_SIZE_AX);
set(groot, 'defaultAxesLineWidth', LINE_WIDTH_AX);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

%% ------------------------- LOAD VALIDATION MAP --------------------------

Sval = load(validationMapFile);
CF = find_first_table(Sval, ["CarryForward_Periods","CarryForward","HSubSupported_PeriodMap"]);

if isempty(CF)
    error('Could not find a carry-forward table in %s.', validationMapFile);
end

CF = standardise_table_types(CF);

% Keep only carry-forward rows if the column exists.
carryCol = find_col(CF, ["CarryForward","UseForAnalysis","HSubSupported","Validated"]);
if ~isempty(carryCol)
    CF = CF(to_logical(CF.(carryCol)), :);
end

% Raw candidates only if available.
srcColCF = find_col(CF, ["Source","SignalSource"]);
if ~isempty(srcColCF)
    CF = CF(strcmpi(string(CF.(srcColCF)), "Raw") | ismissing(string(CF.(srcColCF))), :);
end

% Ultradian bands of interest.
bandColCF = find_col(CF, ["BandName","Band"]);
if isempty(bandColCF)
    error('Validation table does not contain BandName/Band.');
end
CF = CF(ismember(string(CF.(bandColCF)), SELECTED_BANDS), :);

fprintf('Carry-forward rows after filtering selected bands: %d\n', height(CF));

%% ------------------------- LOAD RIDGE PHASE TABLES ----------------------

tsFiles = dir(fullfile(handoffDir, 'WP_TS__*.mat'));
if isempty(tsFiles)
    error('No WP_TS__*.mat files found in selected handoff folder.');
end

allRP = table();

fprintf('Found %d WP_TS files.\n', numel(tsFiles));

for i = 1:numel(tsFiles)
    fpath = fullfile(tsFiles(i).folder, tsFiles(i).name);
    Sts = load(fpath);

    RP = find_first_table(Sts, ["RidgePhase_Long","RidgePhase","RidgePhaseLong"]);

    if isempty(RP)
        warning('No RidgePhase_Long table found in %s. Skipping.', tsFiles(i).name);
        continue;
    end

    RP = standardise_table_types(RP);

    fprintf('Loaded %s (%d RidgePhase rows).\n', tsFiles(i).name, height(RP));

    allRP = [allRP; RP]; %#ok<AGROW>
end

if isempty(allRP)
    error('No RidgePhase_Long rows loaded.');
end

fprintf('\nTotal RidgePhase rows loaded: %d\n', height(allRP));

%% -------------------------- FILTER RIDGE PHASE --------------------------

allRP = standardise_table_types(allRP);

% Required columns in RidgePhase_Long.
requiredRP = ["BandName","RidgePhase_rad"];
for c = requiredRP
    if isempty(find_col(allRP, c))
        error('RidgePhase_Long is missing required column: %s', c);
    end
end

bandColRP = find_col(allRP, ["BandName","Band"]);
sourceColRP = find_col(allRP, ["Source","SignalSource"]);
phaseColRP = find_col(allRP, ["RidgePhase_rad","Phase_rad","WaveletPhase_rad"]);
photoColRP = find_col(allRP, ["Photoperiod_h","PhotoPeriod_h","LightDuration_h"]);
timeColRP = find_col(allRP, ["Time_day","Time_days","Day","TimeDay"]);
ztColRP = find_col(allRP, ["ZT_h","ZT","ZeitgeberTime_h"]);

if isempty(photoColRP)
    error('RidgePhase_Long is missing Photoperiod_h.');
end
if isempty(timeColRP) && isempty(ztColRP)
    error('RidgePhase_Long must contain Time_day or ZT_h.');
end

% Source Raw only.
if ~isempty(sourceColRP)
    allRP = allRP(strcmpi(string(allRP.(sourceColRP)), "Raw"), :);
end

% Selected bands and photoperiods.
allRP = allRP(ismember(string(allRP.(bandColRP)), SELECTED_BANDS), :);
allRP = allRP(ismember(round(double(allRP.(photoColRP)), 6), SELECTED_PHOTOPERIODS), :);

% Validity filters.
if REQUIRE_COI_VALID
    coiCol = find_col(allRP, ["COIValid","COI_valid","ConeValid"]);
    if ~isempty(coiCol)
        allRP = allRP(to_logical(allRP.(coiCol)), :);
    end
end

if REQUIRE_VALID_FLAG
    validCol = find_col(allRP, ["ValidFlag","Valid","PassQC"]);
    if ~isempty(validCol)
        allRP = allRP(to_logical(allRP.(validCol)), :);
    end
end

fprintf('RidgePhase rows after Raw/band/photoperiod/validity filtering: %d\n', height(allRP));

%% ------------------------- JOIN WITH VALIDATION MAP ---------------------

% Prefer CandidateID matching if both tables contain CandidateID.
candColRP = find_col(allRP, ["CandidateID","Candidate_Id","Candidate"]);
candColCF = find_col(CF, ["CandidateID","Candidate_Id","Candidate"]);

fileColRP = find_col(allRP, ["FileStem","File","FileName"]);
fileColCF = find_col(CF, ["FileStem","File","FileName"]);

sigColRP = find_col(allRP, ["SignalID","MouseID","ColumnName","Variable"]);
sigColCF = find_col(CF, ["SignalID","MouseID","ColumnName","Variable"]);

photoColCF = find_col(CF, ["Photoperiod_h","PhotoPeriod_h","LightDuration_h"]);
bandColCF = find_col(CF, ["BandName","Band"]);

if isempty(candColRP) || isempty(candColCF)
    warning(['CandidateID not found in both tables. Falling back to ', ...
        'File/Signal/Photoperiod/Band matching. This is less specific.']);
end

allRP.ValidatedCarryForward = false(height(allRP),1);

if ~isempty(candColRP) && ~isempty(candColCF)
    keyRP = make_key(allRP, {fileColRP, sigColRP, photoColRP, bandColRP, candColRP});
    keyCF = make_key(CF,    {fileColCF, sigColCF, photoColCF, bandColCF, candColCF});
else
    keyRP = make_key(allRP, {fileColRP, sigColRP, photoColRP, bandColRP});
    keyCF = make_key(CF,    {fileColCF, sigColCF, photoColCF, bandColCF});
end

allRP.ValidatedCarryForward = ismember(keyRP, keyCF);
allRP = allRP(allRP.ValidatedCarryForward, :);

fprintf('RidgePhase rows after carry-forward validation filter: %d\n', height(allRP));

if isempty(allRP)
    error('No RidgePhase rows remain after validation filtering. Check CandidateID/key matching.');
end

%% -------------------------- DETECT PHASE EVENTS -------------------------

fprintf('\nDetecting upward 0-radian phase-crossing events...\n');

Events = detect_phase_events(allRP);

if isempty(Events)
    error('No phase events detected. Check phase column and grouping.');
end

fprintf('Detected %d validated ridge-phase events.\n', height(Events));

%% ------------------------- ADD LIGHT STATE LABELS -----------------------

Events.LightState = strings(height(Events),1);
for i = 1:height(Events)
    photo = Events.Photoperiod_h(i);
    zt = Events.EventZT_h(i);

    if photo >= 24
        Events.LightState(i) = "LL";
    elseif zt >= 0 && zt < photo
        Events.LightState(i) = "Light";
    else
        Events.LightState(i) = "Dark";
    end
end

%% --------------------------- SUMMARY TABLES -----------------------------

% Use local grouping helpers rather than groupsummary, because groupsummary
% syntax differs across MATLAB releases/toolbox states.
EventCounts = count_rows_by_group(Events, {'Photoperiod_h','BandName'}, 'N_Events');

CandidateRows = unique(Events(:, {'Photoperiod_h','BandName','CandidateKey'}));
CandidateCounts = count_rows_by_group(CandidateRows, {'Photoperiod_h','BandName'}, ...
    'N_Candidates_WithEvents');

summaryFile = fullfile(outDir, 'ValidatedUR_PhaseEvents_Output.xlsx');
writecell({ ...
    'Validated ultradian ridge-phase event plotting output'; ...
    'Event definition: upward crossing of 0 radians along validated raw ridge-following phase'; ...
    'These are not behavioural activity onsets.'; ...
    ['Generated: ' datestr(now)]}, summaryFile, 'Sheet', 'README', 'Range', 'A1');

writetable(Events, summaryFile, 'Sheet', 'PhaseEvents_Long');
writetable(EventCounts, summaryFile, 'Sheet', 'EventCounts');
writetable(CandidateCounts, summaryFile, 'Sheet', 'CandidateCounts');

%% ----------------------------- ZT HISTOGRAMS ----------------------------

fprintf('\nGenerating ZT histograms...\n');

ztEdges = 0:ZT_BIN_WIDTH_H:24;

for b = 1:numel(SELECTED_BANDS)
    band = SELECTED_BANDS(b);

    for p = 1:numel(SELECTED_PHOTOPERIODS)
        photo = SELECTED_PHOTOPERIODS(p);

        idx = Events.BandName == band & Events.Photoperiod_h == photo;
        if ~any(idx)
            continue;
        end

        E = Events(idx, :);

        f = figure('Color','w','Position',[100 100 1300 700]);
        ax = axes(f); hold(ax, 'on');

        add_light_dark_shading(ax, photo, 0, 24, true, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);

        counts = histcounts(E.EventZT_h, ztEdges);
        centres = ztEdges(1:end-1) + diff(ztEdges)/2;

        bar(ax, centres, counts, 1.0, ...
            'FaceColor', [0.35 0.35 0.35], ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.9);

        add_anchor_lines(ax, photo, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);

        xlim(ax, [0 24]);
        ylim(ax, [0 max(counts)*1.20 + 1]);

        xticks(ax, 0:4:24);
        xlabel(ax, 'ZT / projected ZT (h)', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
        ylabel(ax, 'Validated ridge-phase event count', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');

        title(ax, sprintf('%s | %s | ridge-phase events across ZT', ...
            photo_label(photo), strrep(band, '_', '\_')), ...
            'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');

        style_axes(ax);

        baseName = sprintf('ZT_Hist_%s_%s', photo_label(photo), band);
        save_figure(f, ztHistDir, baseName, FIG_DPI, SAVE_JPEG, SAVE_TIFF);
        close(f);
    end
end

%% ------------------------------ ZT RASTERS ------------------------------

if MAKE_RASTERS
    fprintf('Generating ZT rasters...\n');

    for b = 1:numel(SELECTED_BANDS)
        band = SELECTED_BANDS(b);

        for p = 1:numel(SELECTED_PHOTOPERIODS)
            photo = SELECTED_PHOTOPERIODS(p);

            idx = Events.BandName == band & Events.Photoperiod_h == photo;
            if ~any(idx)
                continue;
            end

            E = Events(idx, :);
            candKeys = unique(E.CandidateKey, 'stable');

            if numel(candKeys) < MIN_EVENTS_PER_GROUP_FOR_RASTER
                continue;
            end

            f = figure('Color','w','Position',[100 100 1400 900]);
            ax = axes(f); hold(ax, 'on');

            add_light_dark_shading(ax, photo, 0, 24, true, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);

            for k = 1:numel(candKeys)
                eidx = E.CandidateKey == candKeys(k);
                zt = E.EventZT_h(eidx);
                y = k * ones(size(zt));

                plot(ax, zt, y, '|', ...
                    'Color', [0 0 0], ...
                    'MarkerSize', 10, ...
                    'LineWidth', 1.2);
            end

            add_anchor_lines(ax, photo, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);

            xlim(ax, [0 24]);
            ylim(ax, [0 numel(candKeys)+1]);
            xticks(ax, 0:4:24);

            xlabel(ax, 'ZT / projected ZT (h)', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
            ylabel(ax, 'Validated ridge candidate', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');

            title(ax, sprintf('%s | %s | validated ridge-phase event raster', ...
                photo_label(photo), strrep(band, '_', '\_')), ...
                'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');

            style_axes(ax);

            baseName = sprintf('ZT_Raster_%s_%s', photo_label(photo), band);
            save_figure(f, rasterDir, baseName, FIG_DPI, SAVE_JPEG, SAVE_TIFF);
            close(f);
        end
    end
end

%% ----------------------- TRANSITION-ALIGNED DENSITIES -------------------

if MAKE_TRANSITION_ALIGNED
    fprintf('Generating transition-aligned event-density plots...\n');

    RelEvents = build_relative_event_table(Events, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);
    writetable(RelEvents, summaryFile, 'Sheet', 'TransitionAligned_Events');

    relEdges = -REL_WINDOW_H_DEFAULT:REL_BIN_WIDTH_H:REL_WINDOW_H_DEFAULT;

    transitionOrder = ["DL","LD","MidLight","MidDark", ...
                       "ProjectedDL_LL","ProjectedLD_LL", ...
                       "ProjectedMidLight_LL","ProjectedMidDark_LL"];

    for b = 1:numel(SELECTED_BANDS)
        band = SELECTED_BANDS(b);

        for p = 1:numel(SELECTED_PHOTOPERIODS)
            photo = SELECTED_PHOTOPERIODS(p);

            idxBase = RelEvents.BandName == band & RelEvents.Photoperiod_h == photo;
            if ~any(idxBase)
                continue;
            end

            transitionsHere = transitionOrder(ismember(transitionOrder, unique(RelEvents.TransitionType(idxBase))));

            for tr = transitionsHere
                idx = idxBase & RelEvents.TransitionType == tr;
                if ~any(idx)
                    continue;
                end

                R = RelEvents(idx, :);
                counts = histcounts(R.RelativeTime_h, relEdges);
                centres = relEdges(1:end-1) + diff(relEdges)/2;

                f = figure('Color','w','Position',[100 100 1200 700]);
                ax = axes(f); hold(ax, 'on');

                bar(ax, centres, counts, 1.0, ...
                    'FaceColor', [0.35 0.35 0.35], ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.9);

                xline(ax, 0, 'k--', 'LineWidth', 1.5);
                xlim(ax, [-REL_WINDOW_H_DEFAULT REL_WINDOW_H_DEFAULT]);
                ylim(ax, [0 max(counts)*1.20 + 1]);

                xlabel(ax, sprintf('Time relative to %s (h)', tr), ...
                    'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
                ylabel(ax, 'Validated ridge-phase event count', ...
                    'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');

                title(ax, sprintf('%s | %s | events aligned to %s', ...
                    photo_label(photo), strrep(band, '_', '\_'), tr), ...
                    'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');

                style_axes(ax);

                baseName = sprintf('RelEvent_%s_%s_%s', photo_label(photo), band, tr);
                save_figure(f, relDir, baseName, FIG_DPI, SAVE_JPEG, SAVE_TIFF);
                close(f);
            end
        end
    end
end

%% ------------------------ COMBINED ZT HISTOGRAM FIGURES -----------------

if MAKE_COMBINED_ZT_FIGURES
    fprintf('Generating combined ZT histogram figures...\n');

    for b = 1:numel(SELECTED_BANDS)
        band = SELECTED_BANDS(b);

        f = figure('Color','w','Position',[100 100 1500 1000]);
        tl = tiledlayout(numel(SELECTED_PHOTOPERIODS), 1, ...
            'TileSpacing','compact', 'Padding','compact');

        maxCount = 0;
        countStore = cell(numel(SELECTED_PHOTOPERIODS),1);

        for p = 1:numel(SELECTED_PHOTOPERIODS)
            photo = SELECTED_PHOTOPERIODS(p);
            idx = Events.BandName == band & Events.Photoperiod_h == photo;
            counts = histcounts(Events.EventZT_h(idx), ztEdges);
            countStore{p} = counts;
            maxCount = max(maxCount, max(counts));
        end

        for p = 1:numel(SELECTED_PHOTOPERIODS)
            photo = SELECTED_PHOTOPERIODS(p);
            ax = nexttile(tl); hold(ax, 'on');

            add_light_dark_shading(ax, photo, 0, 24, true, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);

            counts = countStore{p};
            centres = ztEdges(1:end-1) + diff(ztEdges)/2;

            bar(ax, centres, counts, 1.0, ...
                'FaceColor', [0.35 0.35 0.35], ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.9);

            add_anchor_lines(ax, photo, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H);

            xlim(ax, [0 24]);
            ylim(ax, [0 maxCount*1.20 + 1]);
            xticks(ax, 0:4:24);

            ylabel(ax, photo_label(photo), 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
            style_axes(ax);

            if p == numel(SELECTED_PHOTOPERIODS)
                xlabel(ax, 'ZT / projected ZT (h)', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
            else
                ax.XTickLabel = [];
            end
        end

        title(tl, sprintf('%s | validated ridge-phase events across ZT', strrep(band, '_', '\_')), ...
            'FontName', FONT_NAME, 'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');

        baseName = sprintf('Combined_ZT_Hist_%s', band);
        save_figure(f, combDir, baseName, FIG_DPI, SAVE_JPEG, SAVE_TIFF);
        close(f);
    end
end

%% ------------------------------- FINISH ---------------------------------

fprintf('\nComplete.\n');
fprintf('Output folder:\n%s\n', outDir);
fprintf('Summary workbook:\n%s\n', summaryFile);
diary off;

%% ============================= FUNCTIONS ================================


function tf = is_empty_col(c)
    % Robustly test whether a table-variable-name handle is absent.
    % Handles char '', string "", missing string, and [].
    if isempty(c)
        tf = true;
        return;
    end

    if isstring(c)
        tf = ismissing(c) || strlength(c) == 0;
        return;
    end

    if ischar(c)
        tf = isempty(c);
        return;
    end

    tf = false;
end


function Gtab = count_rows_by_group(T, groupVars, countName)
    % Version-stable replacement for groupsummary(...,'numel',...).
    % Returns one row per unique group and a count column named countName.

    if isempty(T)
        Gtab = table();
        return;
    end

    groupVars = cellstr(groupVars);

    for i = 1:numel(groupVars)
        if ~ismember(groupVars{i}, T.Properties.VariableNames)
            error('Grouping variable missing: %s', groupVars{i});
        end
    end

    G = T(:, groupVars);

    % Convert string/categorical/cellstr group variables to string for stable grouping.
    for i = 1:numel(groupVars)
        v = G.(groupVars{i});
        if iscellstr(v) || isstring(v) || iscategorical(v)
            G.(groupVars{i}) = string(v);
        end
    end

    [Gu, ~, ic] = unique(G, 'rows', 'stable');
    counts = accumarray(ic, 1);

    Gtab = Gu;
    Gtab.(countName) = counts;
end

function make_dir(d)
    if ~exist(d, 'dir')
        mkdir(d);
    end
end

function T = standardise_table_types(T)
    if ~istable(T)
        return;
    end

    vn = T.Properties.VariableNames;

    for i = 1:numel(vn)
        x = T.(vn{i});

        if iscellstr(x) || isstring(x) || iscategorical(x)
            T.(vn{i}) = string(x);
        end
    end
end

function col = find_col(T, candidates)
    col = '';
    if isempty(T) || ~istable(T)
        return;
    end

    vn = string(T.Properties.VariableNames);
    candidates = string(candidates);

    for c = candidates
        idx = strcmpi(vn, c);
        if any(idx)
            col = T.Properties.VariableNames{find(idx,1,'first')};
            return;
        end
    end

    % Try after makeValidName conversion.
    vnValid = string(matlab.lang.makeValidName(T.Properties.VariableNames));
    candValid = string(matlab.lang.makeValidName(cellstr(candidates)));

    for c = candValid
        idx = strcmpi(vnValid, c);
        if any(idx)
            col = T.Properties.VariableNames{find(idx,1,'first')};
            return;
        end
    end
end

function tf = to_logical(x)
    if islogical(x)
        tf = x;
    elseif isnumeric(x)
        tf = x ~= 0 & ~isnan(x);
    else
        xs = lower(strtrim(string(x)));
        tf = xs == "true" | xs == "1" | xs == "yes" | xs == "y" | xs == "pass";
    end

    tf = tf(:);
end

function T = find_first_table(S, preferredNames)
    T = table();

    preferredNames = string(preferredNames);

    % Direct preferred name search.
    for nm = preferredNames
        if isfield(S, nm) && istable(S.(nm))
            T = S.(nm);
            return;
        end
    end

    % Recursive search by preferred name.
    f = fieldnames(S);
    for i = 1:numel(f)
        val = S.(f{i});
        T = recursive_table_search(val, preferredNames);
        if ~isempty(T)
            return;
        end
    end

    % Fallback: first table found.
    for i = 1:numel(f)
        val = S.(f{i});
        T = recursive_first_table(val);
        if ~isempty(T)
            return;
        end
    end
end

function T = recursive_table_search(x, preferredNames)
    T = table();

    if istable(x)
        return;
    end

    if isstruct(x)
        f = fieldnames(x);
        for i = 1:numel(f)
            if any(strcmpi(string(f{i}), preferredNames)) && istable(x.(f{i}))
                T = x.(f{i});
                return;
            end
        end

        for i = 1:numel(f)
            T = recursive_table_search(x.(f{i}), preferredNames);
            if ~isempty(T)
                return;
            end
        end
    end
end

function T = recursive_first_table(x)
    T = table();

    if istable(x)
        T = x;
        return;
    end

    if isstruct(x)
        f = fieldnames(x);
        for i = 1:numel(f)
            T = recursive_first_table(x.(f{i}));
            if ~isempty(T)
                return;
            end
        end
    end
end

function key = make_key(T, cols)
    n = height(T);
    key = strings(n,1);

    for i = 1:numel(cols)
        c = cols{i};

        if is_empty_col(c)
            part = repmat("NA", n, 1);
        else
            x = T.(c);

            if isnumeric(x)
                part = string(round(double(x), 6));
            else
                part = string(x);
            end

            part(ismissing(part)) = "NA";
        end

        if i == 1
            key = part;
        else
            key = key + "|" + part;
        end
    end
end

function Events = detect_phase_events(RP)
    Events = table();

    fileCol = find_col(RP, ["FileStem","File","FileName"]);
    sigCol = find_col(RP, ["SignalID","MouseID","ColumnName","Variable"]);
    photoCol = find_col(RP, ["Photoperiod_h","PhotoPeriod_h","LightDuration_h"]);
    bandCol = find_col(RP, ["BandName","Band"]);
    candCol = find_col(RP, ["CandidateID","Candidate_Id","Candidate"]);
    timeCol = find_col(RP, ["Time_day","Time_days","Day","TimeDay"]);
    ztCol = find_col(RP, ["ZT_h","ZT","ZeitgeberTime_h"]);
    phaseCol = find_col(RP, ["RidgePhase_rad","Phase_rad","WaveletPhase_rad"]);

    if isempty(candCol)
        warning('CandidateID missing. Candidate grouping will use File/Signal/Photoperiod/Band only.');
    end

    groupKey = make_key(RP, {fileCol, sigCol, photoCol, bandCol, candCol});
    RP.GroupKey = groupKey;

    groups = unique(groupKey, 'stable');

    eventRows = {};

    for g = 1:numel(groups)
        idx = RP.GroupKey == groups(g);
        G = RP(idx, :);

        % Sort by time if available.
        if ~isempty(timeCol)
            [~, ord] = sort(double(G.(timeCol)));
        elseif ~isempty(ztCol)
            [~, ord] = sort(double(G.(ztCol)));
        else
            ord = (1:height(G))';
        end
        G = G(ord, :);

        phi = double(G.(phaseCol));
        valid = ~isnan(phi);

        if ~isempty(timeCol)
            tDay = double(G.(timeCol));
        else
            zt = double(G.(ztCol));
            tDay = zt ./ 24;
        end

        valid = valid & ~isnan(tDay);

        phi = phi(valid);
        tDay = tDay(valid);
        Gv = G(valid, :);

        if numel(phi) < 3
            continue;
        end

        phiU = unwrap(phi);
        cyc = floor(phiU ./ (2*pi));
        crossings = find(diff(cyc) > 0);

        if isempty(crossings)
            continue;
        end

        for k = crossings(:)'
            if k >= numel(phiU)
                continue;
            end

            target = 2*pi*cyc(k+1);
            denom = phiU(k+1) - phiU(k);

            if denom == 0 || isnan(denom)
                frac = 0;
            else
                frac = (target - phiU(k)) ./ denom;
                frac = max(0, min(1, frac));
            end

            eDay = tDay(k) + frac*(tDay(k+1) - tDay(k));
            eZT = mod(eDay*24, 24);

            fileVal = get_value_or_na(Gv, fileCol, k);
            sigVal = get_value_or_na(Gv, sigCol, k);
            photoVal = get_value_numeric(Gv, photoCol, k);
            bandVal = string(get_value_or_na(Gv, bandCol, k));
            candVal = get_value_or_na(Gv, candCol, k);

            candKey = string(fileVal) + "|" + string(sigVal) + "|" + ...
                      string(round(photoVal, 6)) + "|" + bandVal + "|" + string(candVal);

            eventRows(end+1, :) = { ...
                string(fileVal), string(sigVal), photoVal, bandVal, string(candVal), ...
                candKey, eDay, eZT, k, "PhaseZeroUpCrossing"}; %#ok<AGROW>
        end
    end

    if isempty(eventRows)
        Events = table();
        return;
    end

    Events = cell2table(eventRows, 'VariableNames', { ...
        'FileStem','SignalID','Photoperiod_h','BandName','CandidateID', ...
        'CandidateKey','EventTime_day','EventZT_h','CrossingIndex','EventType'});
end

function val = get_value_or_na(T, col, row)
    if is_empty_col(col)
        val = "NA";
    else
        x = T.(col);
        if row > numel(x)
            val = "NA";
        else
            val = x(row);
            if iscell(val)
                val = val{1};
            end
        end
    end
end

function val = get_value_numeric(T, col, row)
    if is_empty_col(col)
        val = NaN;
    else
        x = T.(col);
        val = double(x(row));
    end
end

function add_light_dark_shading(ax, photo, xmin, xmax, showProjectedLL, llRefPhoto)
    yl = ylim(ax);

    if photo >= 24
        if showProjectedLL
            % LL has no true dark phase. Add a faint projected former dark
            % window from prior L22 schedule for reference only.
            darkStart = llRefPhoto;
            darkEnd = 24;
            patch(ax, [darkStart darkEnd darkEnd darkStart], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                [0.88 0.88 0.88], ...
                'EdgeColor','none', ...
                'FaceAlpha',0.35);
        end
    else
        darkStart = photo;
        darkEnd = 24;
        patch(ax, [darkStart darkEnd darkEnd darkStart], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.85 0.85], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.45);
    end

    uistack(findobj(ax,'Type','patch'), 'bottom');
    xlim(ax, [xmin xmax]);
end

function add_anchor_lines(ax, photo, llRefPhoto)
    if photo >= 24
        xline(ax, 0, '--', 'Projected DL', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, ...
            'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal');
        xline(ax, llRefPhoto, '--', 'Projected LD', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, ...
            'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal');
    else
        xline(ax, 0, '--', 'DL', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, ...
            'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal');
        xline(ax, photo, '--', 'LD', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, ...
            'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal');
    end
end

function RelEvents = build_relative_event_table(Events, llRefPhoto)
    rows = {};

    for i = 1:height(Events)
        photo = Events.Photoperiod_h(i);
        zt = Events.EventZT_h(i);

        if photo >= 24
            anchors = table( ...
                ["ProjectedDL_LL"; "ProjectedLD_LL"; "ProjectedMidLight_LL"; "ProjectedMidDark_LL"], ...
                [0; llRefPhoto; llRefPhoto/2; llRefPhoto + ((24-llRefPhoto)/2)], ...
                'VariableNames', {'TransitionType','AnchorZT_h'});
        else
            anchors = table( ...
                ["DL"; "LD"; "MidLight"; "MidDark"], ...
                [0; photo; photo/2; photo + ((24-photo)/2)], ...
                'VariableNames', {'TransitionType','AnchorZT_h'});
        end

        for a = 1:height(anchors)
            rel = circular_hour_diff(zt, anchors.AnchorZT_h(a));

            rows(end+1,:) = { ...
                Events.FileStem(i), Events.SignalID(i), photo, Events.BandName(i), ...
                Events.CandidateID(i), Events.CandidateKey(i), Events.EventTime_day(i), ...
                zt, anchors.TransitionType(a), anchors.AnchorZT_h(a), rel}; %#ok<AGROW>
        end
    end

    RelEvents = cell2table(rows, 'VariableNames', { ...
        'FileStem','SignalID','Photoperiod_h','BandName','CandidateID', ...
        'CandidateKey','EventTime_day','EventZT_h','TransitionType', ...
        'AnchorZT_h','RelativeTime_h'});
end

function d = circular_hour_diff(t, anchor)
    d = mod(t - anchor + 12, 24) - 12;
end

function txt = photo_label(photo)
    if photo >= 24
        txt = "LL";
    else
        txt = "L" + string(round(photo));
    end
end

function style_axes(ax)
    set(ax, 'FontName', 'Times New Roman', ...
        'FontSize', 18, ...
        'LineWidth', 1.5, ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'Layer', 'top');
end

function save_figure(f, outDir, baseName, dpi, saveJpeg, saveTiff)
    if saveJpeg
        jpgFile = fullfile(outDir, [char(baseName) '.jpg']);
        exportgraphics(f, jpgFile, 'Resolution', dpi);
        fprintf('Saved: %s\n', jpgFile);
    end

    if saveTiff
        tifFile = fullfile(outDir, [char(baseName) '.tif']);
        exportgraphics(f, tifFile, 'Resolution', dpi);
        fprintf('Saved: %s\n', tifFile);
    end
end