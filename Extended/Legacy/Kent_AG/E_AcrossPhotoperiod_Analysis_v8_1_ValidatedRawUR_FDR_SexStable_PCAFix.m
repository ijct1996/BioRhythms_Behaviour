%% ========================================================================
% Script 5 (VALIDATED RAW UR + LME BH/FDR + optional sex/stable-day safeguards, v8.0)
% File: AcrossPhotoperiod_Analysis_v8_ValidatedRawUR_FDR_SexStable.m
% ========================================================================
% Purpose
%   Load handoff mats from Script 1 (e.g., Behav_wavelet_v9.m):
%     - WP_Summary__*.mat (fast, summary-first)
%     - WP_TS__*.mat (slow, streamed last)
%     - Optional WP_Package_Index.xlsx (preferred TS mapping)
%
% v7 conceptual change
%   - Final biological characterisation uses Raw data for both circadian and
%     ultradian bands.
%   - Selective-HSub is no longer used as the ultradian signal source. Instead,
%     the HSubSupported_PeriodMap from Script 3 is used to retain only Raw
%     ultradian candidates that survived Selective-HSub validation.
%   - HSub/Full-Ladder outputs are retained as validation and sensitivity metadata.
%
% v7.1 inferential update
%   - Mixed-effects coefficient and ANOVA p-values are collated into a
%     transparent LME inference table.
%   - Benjamini-Hochberg FDR correction is applied within predefined
%     families: Delta_LME_Coefficients, Delta_LME_Anova,
%     Power_LME_Coefficients and Power_LME_Anova.
%   - Intercepts are retained in raw coefficient outputs but excluded from
%     BH/FDR correction because they are not biological hypothesis tests.
%
% Key features in v6.3 (carried into v6.3.1)
%   1) Colour governance (publication-grade consistency):
%       - Paul Tol "Bright" palette is the canonical palette.
%       - Band colours: each band (CR + each UR) gets a unique colour, used
%         consistently across ALL band-encoded outputs (CRvsUR, composition).
%       - Photoperiod colours: each photoperiod condition gets a unique colour
%         (Tol Bright). This mapping is used consistently across ALL
%         photoperiod-encoded bar outputs and PCA.
%       - If Script 1 handoff provides valid colour maps, they are used; otherwise
%         robust defaults are assigned internally.
%
%   2) Global plotting consistency:
%       - No top x-axis and no right y-axis (axes box off)
%       - Tick marks outside
%       - Times New Roman everywhere
%       - Titles and axis labels bold
%       - JPEG only, 600 DPI only
%
%   3) LME handling:
%       - LME outputs are written to a SEPARATE workbook only if at least one
%         model successfully fits AND produces table outputs.
%       - If all LME work fails: no LME folders/files, log only.
%
%   4) Whole-experiment scalograms:
%       - Removed entirely (folders + code + settings)
%
% Fix in v6.3.1
%   - PCA indexing bugfix:
%       Old (wrong): ppU = unique(pp); ppU = sort(ppU(~isnan(pp)));
%       New (right): ppU = unique(pp); ppU = sort(ppU(~isnan(ppU)));
%
% MATLAB: R2025a/R2025b
% ========================================================================

clearvars; close all; clc;

%% ----------------------------- SETTINGS ---------------------------------
SAVE_DPI  = 600;
FONT_NAME = 'Times New Roman';
TICK_DIR  = 'out';

PHASE_ORDER = {'All','Light','Dark'};

CR_BAND  = "CR_20_28";
UR_BANDS = ["UR_1_3","UR_3_6","UR_6_9","UR_9_12","UR_12_18"];
BANDS_ALL = [CR_BAND UR_BANDS];

SRC_CR = "Raw";
SRC_UR = "Raw";

% v7 validated-Raw ultradian settings
USE_HSUB_VALIDATION_MAP = true;
REQUIRE_VALIDATION_MAP  = true;
VALIDATION_PRIMARY_PHASE = "All";
VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK = true;

% Robust plotting defaults
USE_ROBUST_RATIO = true;  % median + IQR instead of mean ± SD for UR/CR ratios

% Summary-level inference
DO_LME = true;
LME_USE_PHASE = true;

% Multiple-testing correction for inferential outputs
DO_LME_BH_FDR = true;
FDR_ALPHA = 0.05;

% v8 sex and stable-day safeguards
INCLUDE_SEX_IN_LME_IF_AVAILABLE = true;   % Adds Sex as an additive fixed effect if a usable Sex column can be found or inferred
ALLOW_SEX_INFERENCE_FROM_SIGNALID = true; % Cautiously infer M/F only from clear SignalID labels if no Sex column exists
SEX_FIXED_EFFECT_MODE = "AdditiveOnly";   % Avoids high-order interactions unless explicitly introduced later

% Stable-day filter for streamed time-series analyses. Summary-level tables from
% Behav_wavelet are already aggregated, so this filter applies to streamed TS
% outputs here. For summary-level stable-day inference, regenerate Script 2
% summaries with an equivalent stable-day option.
EXCLUDE_INITIAL_DAYS_FROM_TS_ANALYSES = false;
MIN_DAY_FOR_STABLE_TS = 0.0;

% PCA settings
DO_PCA_POWERS_ONLY = true;   % CR + UR band powers (avoid redundancy)
DO_PCA_DELTAS_ONLY = true;   % optional second PCA on deltas only
DO_PCA_CLUSTERING  = true;   % produce cluster-highlight figs + assignment tables

% Time-series analyses (streamed; last)
DO_DOMINANCE_BP   = true;    % P(UR_b > CR) using BandPower_Long
DO_COVERAGE       = true;    % coverage fraction from ValidFlag
DO_STABILITY_RP   = true;    % stability of RidgePeriod over time (SD + MAD)
DO_STABILITY_RPOW = true;    % stability of RidgePower over time (SD + MAD)

TS_PHASES = {'All'};         % keep robust; extend later

NEED_TS = DO_DOMINANCE_BP || DO_COVERAGE || DO_STABILITY_RP || DO_STABILITY_RPOW;

%% ---------------------- SELECT HANDOFF ROOT -----------------------------
rootSel = uigetdir(pwd, 'Select folder containing WP_Summary__*.mat (handoff folder or its parent)');
if isequal(rootSel,0)
    fprintf('No folder selected. Exiting.\n');
    return;
end

ddS = dir(fullfile(rootSel, '**', 'WP_Summary__*.mat'));
if isempty(ddS)
    fprintf('No WP_Summary__*.mat found under:\n  %s\n', rootSel);
    return;
end

sumPaths = cell(numel(ddS),1);
for i = 1:numel(ddS)
    sumPaths{i} = fullfile(ddS(i).folder, ddS(i).name);
end

fprintf('Searching for summary packages under:\n  %s\n', rootSel);
fprintf('Discovered %d summary package(s).\n', numel(sumPaths));

%% ---------------------- SELECT OUTPUT FOLDER ----------------------------
outDir = uigetdir(pwd, 'Select output folder for AcrossPhotoperiod results');
if isequal(outDir,0)
    fprintf('No output folder selected. Exiting.\n');
    return;
end

dirs = struct();
dirs.Run    = outDir;
dirs.Tables = fullfile(outDir, 'Tables');
dirs.Fig    = fullfile(outDir, 'Figures');
dirs.Logs   = fullfile(outDir, 'Logs');

dirs.Fig_CRonly   = fullfile(dirs.Fig, 'CR_Absolute');
dirs.Fig_URonly   = fullfile(dirs.Fig, 'UR_Absolute');
dirs.Fig_CRvsUR   = fullfile(dirs.Fig, 'CR_vs_UR_Bandwise');
dirs.Fig_Delta    = fullfile(dirs.Fig, 'DeltaLog10_URminusCR');
dirs.Fig_Ratio    = fullfile(dirs.Fig, 'Ratio_URoverCR');
dirs.Fig_URComp   = fullfile(dirs.Fig, 'UR_Composition_FracLinear');
dirs.Fig_CompAll  = fullfile(dirs.Fig, 'Composition_CRplusUR_FracLinear');
dirs.Fig_PCA      = fullfile(dirs.Fig, 'PCA');

dirs.Fig_Dom      = fullfile(dirs.Fig, 'Dominance_Occupancy_BP');
dirs.Fig_Cov      = fullfile(dirs.Fig, 'Coverage');
dirs.Fig_StabRP   = fullfile(dirs.Fig, 'Stability_RidgePeriod');
dirs.Fig_StabRPOW = fullfile(dirs.Fig, 'Stability_RidgePower');

ensure_dir(dirs.Run); ensure_dir(dirs.Tables); ensure_dir(dirs.Fig); ensure_dir(dirs.Logs);
ensure_dir(dirs.Fig_CRonly); ensure_dir(dirs.Fig_URonly); ensure_dir(dirs.Fig_CRvsUR);
ensure_dir(dirs.Fig_Delta); ensure_dir(dirs.Fig_Ratio); ensure_dir(dirs.Fig_URComp);
ensure_dir(dirs.Fig_CompAll);
ensure_dir(dirs.Fig_PCA);
ensure_dir(dirs.Fig_Dom); ensure_dir(dirs.Fig_Cov); ensure_dir(dirs.Fig_StabRP); ensure_dir(dirs.Fig_StabRPOW);

xlsxOut = fullfile(dirs.Tables, 'AcrossPhotoperiod_Outputs.xlsx');
if exist(xlsxOut,'file'), delete(xlsxOut); end

runLog = fullfile(dirs.Logs, 'Run_Log.txt');
write_text(runLog, sprintf('AcrossPhotoperiod run started: %s\nRoot: %s\nOut: %s\n', datestr(now,31), rootSel, outDir));

writecell({'AcrossPhotoperiod outputs. Summary sheets written first; TS sheets appended later.'}, xlsxOut, 'Sheet', 'README');

%% ---------------- LOAD HSUB-SUPPORTED PERIOD MAP -----------------------
% The v7 pipeline requires the Script 3 carry-forward map. This map defines
% which Raw ultradian candidates survived Selective-HSub validation. The map
% is used as a gate for all final ultradian analyses below.
ValidationMap = struct();
CarryForward_Periods = table();
validationMapPath = '';

if USE_HSUB_VALIDATION_MAP
    % Try to auto-locate the map first.
    mapHits = dir(fullfile(rootSel, '**', 'HSubSupported_PeriodMap.mat'));
    if ~isempty(mapHits)
        validationMapPath = fullfile(mapHits(1).folder, mapHits(1).name);
    else
        [vf, vp] = uigetfile({'HSubSupported_PeriodMap.mat','HSubSupported_PeriodMap.mat'}, ...
            'Select HSubSupported_PeriodMap.mat from Script 3');
        if ~isequal(vf,0)
            validationMapPath = fullfile(vp, vf);
        end
    end

    if isempty(validationMapPath) || ~isfile(validationMapPath)
        if REQUIRE_VALIDATION_MAP
            fprintf('No HSubSupported_PeriodMap.mat selected/found. Exiting because REQUIRE_VALIDATION_MAP=true.\n');
            return;
        else
            warning('No validation map found. Continuing as unvalidated Raw-UR analysis.');
        end
    else
        VM = load(validationMapPath);
        if isfield(VM,'validationMap')
            ValidationMap = VM.validationMap;
        else
            ValidationMap = VM;
        end

        if isfield(ValidationMap,'CarryForward_Periods')
            CarryForward_Periods = ValidationMap.CarryForward_Periods;
        elseif isfield(VM,'CarryForward_Periods')
            CarryForward_Periods = VM.CarryForward_Periods;
        else
            error('Validation map loaded, but no CarryForward_Periods table was found.');
        end

        CarryForward_Periods = harmonise_carryforward_table(CarryForward_Periods);
        fprintf('Loaded validation map: %s\n', validationMapPath);
        fprintf('Carry-forward validated Raw UR candidates: %d\n', height(CarryForward_Periods));
    end
end

append_text(runLog, sprintf('Validation map path: %s\nCarry-forward rows: %d\n', validationMapPath, height(CarryForward_Periods)));

% LME log only (always available)
lmeLog = fullfile(dirs.Logs, 'LME_Log.txt');
delete_if_exists(lmeLog);
LOG = fopen(lmeLog, 'a');
fprintf(LOG, 'LME log opened: %s\n', datestr(now,31));

%% ========================================================================
% Resolve WP_Package_Index.xlsx (optional but preferred)
% ========================================================================
IndexTable = table();
idxHits = dir(fullfile(rootSel, '**', 'WP_Package_Index.xlsx'));
if ~isempty(idxHits)
    idxPath = fullfile(idxHits(1).folder, idxHits(1).name);
    try
        IndexTable = readtable(idxPath, 'Sheet', 'Index', 'VariableNamingRule','preserve');
        fprintf('Found index: %s\n', idxPath);
        append_text(runLog, sprintf('Found WP_Package_Index.xlsx: %s\n', idxPath));
    catch ME
        IndexTable = table();
        append_text(runLog, sprintf('Index load failed: %s\n', ME.message));
    end
else
    fprintf('No WP_Package_Index.xlsx found under root. TS paths will be resolved by filename search.\n');
    append_text(runLog, 'No WP_Package_Index.xlsx found. Using fallback TS search by filename.\n');
end

%% ========================================================================
% PHASE A: SUMMARY-ONLY
% ========================================================================
fprintf('\nLoading summary packages...\n');

BCS = table();
SumMeta = table();

TSMap = containers.Map('KeyType','char','ValueType','char');
SumMap = containers.Map('KeyType','char','ValueType','char');

% Optional colour payload from Script 1 (if present)
HandoffBandColour = [];
HandoffPPColour   = [];

for i = 1:numel(sumPaths)
    p = sumPaths{i};
    try
        S = load(p, 'pkgS');
        if ~isfield(S,'pkgS') || isempty(S.pkgS), continue; end
        pkgS = S.pkgS;

        % Capture handoff colour maps if present (first hit wins)
        if isempty(HandoffBandColour) || isempty(HandoffPPColour)
            [bc, pc] = try_get_colour_maps_from_pkgS(pkgS);
            if isempty(HandoffBandColour) && ~isempty(bc), HandoffBandColour = bc; end
            if isempty(HandoffPPColour)   && ~isempty(pc), HandoffPPColour   = pc; end
        end

        if isfield(pkgS,'tables') && isfield(pkgS.tables,'BandConditionSummary') && ~isempty(pkgS.tables.BandConditionSummary)
            T = pkgS.tables.BandConditionSummary;
            T.PackageID  = repmat(i, height(T), 1);
            T.SummaryMat = repmat(string(p), height(T), 1);

            if isfield(pkgS,'file') && isfield(pkgS.file,'FileStem')
                T.FileStemPkg = repmat(string(pkgS.file.FileStem), height(T), 1);
            else
                T.FileStemPkg = repmat(missing, height(T), 1);
            end

            BCS = [BCS; T]; %#ok<AGROW>
        end

        stem = extract_stem_from_summary(p);
        if ~isempty(stem)
            if ~isKey(SumMap, stem)
                SumMap(stem) = char(string(p));
            end
        end

        % Resolve TS mat path for this stem
        if NEED_TS && ~isempty(stem)
            tsPath = '';

            % 1) Prefer IndexTable if available
            if ~isempty(IndexTable) && ismember('FileStem', IndexTable.Properties.VariableNames) && ismember('TSMat', IndexTable.Properties.VariableNames)
                ix = find(strcmp(string(IndexTable.FileStem), string(stem)), 1, 'first');
                if ~isempty(ix)
                    tsPath = char(string(IndexTable.TSMat(ix)));
                end
            end

            % 2) Fallback: search under rootSel for WP_TS__<stem>.mat
            if isempty(tsPath) || ~isfile(tsPath)
                ddTS = dir(fullfile(rootSel, '**', ['WP_TS__' stem '.mat']));
                if ~isempty(ddTS)
                    tsPath = fullfile(ddTS(1).folder, ddTS(1).name);
                end
            end

            if ~isempty(tsPath) && isfile(tsPath)
                if ~isKey(TSMap, stem)
                    TSMap(stem) = tsPath;
                end
            else
                append_text(runLog, sprintf('TS not resolved for stem %s (summary: %s)\n', stem, p));
            end
        end

        % Meta
        try
            m = struct();
            if isfield(pkgS,'meta'), m = pkgS.meta; end
            r = table(string(p), string(getfield_or(m,'RunType','')), string(getfield_or(m,'HSubMode','')), ...
                'VariableNames', {'SummaryMat','RunType','HSubMode'});
            SumMeta = [SumMeta; r]; %#ok<AGROW>
        catch
        end

    catch ME
        append_text(runLog, sprintf('Summary load failed: %s\n', ME.message));
    end
end

if isempty(BCS)
    fprintf('No BandConditionSummary found. Exiting.\n');
    fclose(LOG);
    return;
end

fprintf('Summary load complete. BCS rows: %d\n', height(BCS));

BCS.File     = string(BCS.File);
BCS.SignalID = string(BCS.SignalID);
BCS.Source   = string(BCS.Source);
BCS.Phase    = string(BCS.Phase);
BCS.BandName = string(BCS.BandName);
if ismember('FileStemPkg', BCS.Properties.VariableNames)
    BCS.FileStemPkg = string(BCS.FileStemPkg);
else
    BCS.FileStemPkg = repmat(missing, height(BCS), 1);
end

% v7: annotate and filter summary rows for validated-Raw ultradian analyses.
BCS_AllRows = BCS;
if USE_HSUB_VALIDATION_MAP && ~isempty(CarryForward_Periods)
    BCS_AllRows = annotate_BCS_with_validation(BCS_AllRows, CarryForward_Periods, UR_BANDS, SRC_UR, VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK);
    BCS = filter_BCS_for_validated_raw_analysis(BCS_AllRows, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, true);
else
    BCS_AllRows = annotate_BCS_with_validation(BCS_AllRows, table(), UR_BANDS, SRC_UR, VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK);
    BCS = filter_BCS_for_validated_raw_analysis(BCS_AllRows, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, false);
end

fprintf('Summary rows after v7 validated-Raw filter: %d\n', height(BCS));
append_text(runLog, sprintf('Summary rows all: %d | analysis-used rows: %d\n', height(BCS_AllRows), height(BCS)));

if isempty(BCS)
    fprintf('No rows available after validated-Raw filtering. Check CarryForward_Periods and BCS keys. Exiting.\n');
    fclose(LOG);
    return;
end


phasesPresent = unique(BCS.Phase);
phaseList = intersect(string(PHASE_ORDER), phasesPresent, 'stable');
if isempty(phaseList), phaseList = phasesPresent; end

% Photoperiod list (global) and colour maps
ppValsAll = unique(double(BCS.Photoperiod_h));
ppValsAll = sort(ppValsAll(isfinite(ppValsAll)));

% Build canonical colour maps (handoff preferred; fallback to internal Tol Bright)
BandColourMap = build_band_colour_map(BANDS_ALL, HandoffBandColour);
PhotoperiodColourMap = build_photoperiod_colour_map(ppValsAll, HandoffPPColour);

% Write summary sheets early
V7Settings = table( ...
    ["SRC_CR"; "SRC_UR"; "USE_HSUB_VALIDATION_MAP"; "REQUIRE_VALIDATION_MAP"; "VALIDATION_PRIMARY_PHASE"; "VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK"; "DO_LME_BH_FDR"; "FDR_ALPHA"; "INCLUDE_SEX_IN_LME_IF_AVAILABLE"; "ALLOW_SEX_INFERENCE_FROM_SIGNALID"; "SEX_FIXED_EFFECT_MODE"; "EXCLUDE_INITIAL_DAYS_FROM_TS_ANALYSES"; "MIN_DAY_FOR_STABLE_TS"; "ValidationMapPath"], ...
    [string(SRC_CR); string(SRC_UR); string(USE_HSUB_VALIDATION_MAP); string(REQUIRE_VALIDATION_MAP); string(VALIDATION_PRIMARY_PHASE); string(VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK); string(DO_LME_BH_FDR); string(FDR_ALPHA); string(INCLUDE_SEX_IN_LME_IF_AVAILABLE); string(ALLOW_SEX_INFERENCE_FROM_SIGNALID); string(SEX_FIXED_EFFECT_MODE); string(EXCLUDE_INITIAL_DAYS_FROM_TS_ANALYSES); string(MIN_DAY_FOR_STABLE_TS); string(validationMapPath)], ...
    'VariableNames', {'Setting','Value'});
writetable(V7Settings, xlsxOut, 'Sheet', 'Analysis_Settings');
writetable(BCS_AllRows, xlsxOut, 'Sheet', 'BandConditionSummary_All');
writetable(BCS, xlsxOut, 'Sheet', 'BandCondSummary_AnalysisUsed');
if ~isempty(CarryForward_Periods), writetable(CarryForward_Periods, xlsxOut, 'Sheet', 'ValidationMap_Used'); end
if ismember('ValidatedRawUR', BCS_AllRows.Properties.VariableNames)
    RawUR_NotValidated = BCS_AllRows(BCS_AllRows.Source==SRC_UR & ismember(BCS_AllRows.BandName, UR_BANDS) & ~BCS_AllRows.ValidatedRawUR, :);
    if ~isempty(RawUR_NotValidated), writetable(RawUR_NotValidated, xlsxOut, 'Sheet', 'UR_Raw_NotValidated'); end
end
if ~isempty(SumMeta), writetable(SumMeta, xlsxOut, 'Sheet', 'SummaryMeta'); end

%% ---------------- SUMMARY DERIVATIONS -----------------------------------
[Tpair, TpairSummary] = build_CR_UR_pairs_from_BCS(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
Tabs   = build_absolute_summary_tables(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
TURcomp= build_UR_composition_linearfrac(BCS, UR_BANDS, SRC_UR);

% CR+UR composition
TCompAll = build_CRplusUR_composition_linearfrac(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);

writetable(Tpair, xlsxOut, 'Sheet', 'CR_UR_Pairs_PerMouse');
writetable(TpairSummary, xlsxOut, 'Sheet', 'CR_UR_Pairs_Summary');
writetable(Tabs, xlsxOut, 'Sheet', 'AbsolutePower_Summary');
writetable(TURcomp, xlsxOut, 'Sheet', 'UR_Composition_FracLinear');
if ~isempty(TCompAll)
    writetable(TCompAll, xlsxOut, 'Sheet', 'Composition_CRplusUR_FracLinear');
end

fprintf('Wrote summary tables to:\n  %s\n', xlsxOut);

%% ---------------- FIGURES: SUMMARY --------------------------------------
% Absolute CR: photoperiod-encoded colours
plot_absolute_band_across_pp( ...
    Tabs(Tabs.BandName==CR_BAND,:), ...
    dirs.Fig_CRonly, SAVE_DPI, FONT_NAME, TICK_DIR, ...
    'Circadian band (CR_20_28) power across photoperiod', ...
    PhotoperiodColourMap);

% Absolute URs: photoperiod-encoded colours
for b = 1:numel(UR_BANDS)
    ub = UR_BANDS(b);
    plot_absolute_band_across_pp( ...
        Tabs(Tabs.BandName==ub,:), ...
        dirs.Fig_URonly, SAVE_DPI, FONT_NAME, TICK_DIR, ...
        sprintf('Validated Raw ultradian band (%s) power across photoperiod', char(ub)), ...
        PhotoperiodColourMap);
end

% CR vs UR (band-encoded colours), Delta/Ratio (photoperiod-encoded colours)
for ph = phaseList(:).'
    phaseTag = char(ph);
    Dpair = Tpair(Tpair.Phase==string(phaseTag), :);
    if isempty(Dpair), continue; end

    for b = 1:numel(UR_BANDS)
        ub = UR_BANDS(b);
        Db = Dpair(Dpair.UR_Band==ub, :);
        if isempty(Db), continue; end

        S = summarise_CR_UR(Db);
        plot_CR_vs_UR(S, CR_BAND, ub, phaseTag, dirs.Fig_CRvsUR, SAVE_DPI, FONT_NAME, TICK_DIR, BandColourMap);

        Sd = summarise_delta(Db);
        plot_delta(Sd, ub, phaseTag, dirs.Fig_Delta, SAVE_DPI, FONT_NAME, TICK_DIR, PhotoperiodColourMap);

        Sr = summarise_ratio(Db, USE_ROBUST_RATIO);
        plot_ratio(Sr, ub, phaseTag, dirs.Fig_Ratio, SAVE_DPI, FONT_NAME, TICK_DIR, USE_ROBUST_RATIO, PhotoperiodColourMap);
    end
end

% UR-only composition plot: band colours
plot_UR_composition(TURcomp, UR_BANDS, phaseList, dirs.Fig_URComp, SAVE_DPI, FONT_NAME, TICK_DIR, BandColourMap);

% CR+UR composition plot: band colours
if ~isempty(TCompAll)
    plot_CRplusUR_composition(TCompAll, CR_BAND, UR_BANDS, phaseList, dirs.Fig_CompAll, SAVE_DPI, FONT_NAME, TICK_DIR, BandColourMap);
end

%% ---------------- LME (hardened, separate xlsx only if success) ----------
LME_Tables = struct();
LME_Inference_Raw = table();
LME_SuccessCount = 0;

if DO_LME
    fprintf('\nRunning mixed-effects models (summary-level)...\n');
    fprintf(LOG, 'Starting LME run at %s\n', datestr(now,31));

    try
        TpairL = Tpair;
        assert(istable(TpairL), 'Tpair is not a table');

        TpairL.SignalID = categorical(TpairL.SignalID);
        TpairL.File     = categorical(TpairL.File);
        TpairL.Phase    = categorical(TpairL.Phase);
        TpairL.UR_Band  = categorical(TpairL.UR_Band);
        [TpairL, SexInfo_Delta] = add_sex_column_if_available(TpairL, INCLUDE_SEX_IN_LME_IF_AVAILABLE, ALLOW_SEX_INFERENCE_FROM_SIGNALID);
        LME_Tables.SexInfo_Delta = SexInfo_Delta;

        % 1) Delta models (per UR band)
        for b = 1:numel(UR_BANDS)
            ub = UR_BANDS(b);
            D = TpairL(TpairL.UR_Band==categorical(ub), :);
            if has_usable_sex(D)
                D = D(string(D.Sex) ~= "Unknown" & ~ismissing(string(D.Sex)), :);
            end
            if ~istable(D) || height(D) < 3
                fprintf(LOG, 'Delta LME skipped %s: not enough rows\n', char(ub));
                continue;
            end

            hasPhase = LME_USE_PHASE && numel(categories(D.Phase)) > 1;
            hasSex = has_usable_sex(D);
            fml = build_lme_formula('Delta_log10', hasPhase, hasSex, SEX_FIXED_EFFECT_MODE);

            try
                mdl = fitlme(D, fml);
                Tcoef = mdl.Coefficients;
                Tanova = anova(mdl);
                LME_Tables.(safe_field("Coef_Delta_" + ub))  = Tcoef;
                LME_Tables.(safe_field("Anova_Delta_" + ub)) = Tanova;
                LME_Inference_Raw = append_table_rows(LME_Inference_Raw, annotate_lme_inference_table(Tcoef,  "Delta_LME_Coefficients", "Delta", ub, "Delta_log10", fml, "Coefficients"));
                LME_Inference_Raw = append_table_rows(LME_Inference_Raw, annotate_lme_inference_table(Tanova, "Delta_LME_Anova",        "Delta", ub, "Delta_log10", fml, "ANOVA"));
                LME_SuccessCount = LME_SuccessCount + 1;
            catch ME
                fprintf(LOG, 'Delta LME failed %s: %s\n', char(ub), ME.message);
            end
        end

        % 2) Power models (CR + UR bands)
        TBP = build_bandpower_table_from_BCS(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
        assert(istable(TBP), 'TBP is not a table');

        TBP.SignalID = categorical(TBP.SignalID);
        TBP.File     = categorical(TBP.File);
        TBP.Phase    = categorical(TBP.Phase);
        TBP.BandName = categorical(TBP.BandName);
        [TBP, SexInfo_Power] = add_sex_column_if_available(TBP, INCLUDE_SEX_IN_LME_IF_AVAILABLE, ALLOW_SEX_INFERENCE_FROM_SIGNALID);
        LME_Tables.SexInfo_Power = SexInfo_Power;

        for b = 1:numel(BANDS_ALL)
            bn = BANDS_ALL(b);
            D = TBP(TBP.BandName==categorical(bn), :);
            if has_usable_sex(D)
                D = D(string(D.Sex) ~= "Unknown" & ~ismissing(string(D.Sex)), :);
            end
            if ~istable(D) || height(D) < 3
                fprintf(LOG, 'Power LME skipped %s: not enough rows\n', char(bn));
                continue;
            end

            hasPhase = LME_USE_PHASE && numel(categories(D.Phase)) > 1;
            hasSex = has_usable_sex(D);
            fml = build_lme_formula('MeanBandPower_log10', hasPhase, hasSex, SEX_FIXED_EFFECT_MODE);

            try
                mdl = fitlme(D, fml);
                Tcoef = mdl.Coefficients;
                Tanova = anova(mdl);
                LME_Tables.(safe_field("Coef_Power_" + bn))  = Tcoef;
                LME_Tables.(safe_field("Anova_Power_" + bn)) = Tanova;
                LME_Inference_Raw = append_table_rows(LME_Inference_Raw, annotate_lme_inference_table(Tcoef,  "Power_LME_Coefficients", "Power", bn, "MeanBandPower_log10", fml, "Coefficients"));
                LME_Inference_Raw = append_table_rows(LME_Inference_Raw, annotate_lme_inference_table(Tanova, "Power_LME_Anova",        "Power", bn, "MeanBandPower_log10", fml, "ANOVA"));
                LME_SuccessCount = LME_SuccessCount + 1;
            catch ME
                fprintf(LOG, 'Power LME failed %s: %s\n', char(bn), ME.message);
            end
        end

    catch ME
        fprintf(LOG, 'Top-level LME block failed: %s\n', ME.message);
    end

    % Write LME only if any success and table outputs exist
    if LME_SuccessCount > 0 && ~isempty(fieldnames(LME_Tables))
        lmeXlsx = fullfile(dirs.Tables, 'AcrossPhotoperiod_LME_Outputs.xlsx');
        delete_if_exists(lmeXlsx);

        try
            fns = fieldnames(LME_Tables);
            for i2 = 1:numel(fns)
                nm = fns{i2};
                T = LME_Tables.(nm);
                if istable(T)
                    writetable(T, lmeXlsx, 'Sheet', safe_sheet_name(nm));
                else
                    fprintf(LOG, 'Skipping non-table LME output %s\n', nm);
                end
            end

            % v7.1: consolidated inferential table with BH/FDR correction.
            if DO_LME_BH_FDR && istable(LME_Inference_Raw) && ~isempty(LME_Inference_Raw)
                LME_Inference_FDR = apply_bh_fdr_by_family(LME_Inference_Raw, FDR_ALPHA);
                writetable(LME_Inference_Raw, lmeXlsx, 'Sheet', 'LME_Inference_Raw');
                writetable(LME_Inference_FDR, lmeXlsx, 'Sheet', 'LME_Inference_BH_FDR');

                write_lme_fdr_subset(LME_Inference_FDR, lmeXlsx, "Delta", "ANOVA",        'LME_Anova_Delta_BH_FDR');
                write_lme_fdr_subset(LME_Inference_FDR, lmeXlsx, "Power", "ANOVA",        'LME_Anova_Power_BH_FDR');
                write_lme_fdr_subset(LME_Inference_FDR, lmeXlsx, "Delta", "Coefficients", 'LME_Coef_Delta_BH_FDR');
                write_lme_fdr_subset(LME_Inference_FDR, lmeXlsx, "Power", "Coefficients", 'LME_Coef_Power_BH_FDR');
                if ismember('UR12_18SensitivityFlag', LME_Inference_FDR.Properties.VariableNames)
                    T12 = LME_Inference_FDR(logical(LME_Inference_FDR.UR12_18SensitivityFlag), :);
                    if ~isempty(T12), writetable(T12, lmeXlsx, 'Sheet', 'UR12_18_Sensitivity_FDR'); end
                end
            end

            fprintf(LOG, 'LME workbook written: %s\n', lmeXlsx);
        catch ME
            fprintf(LOG, 'Writing LME workbook failed: %s\n', ME.message);
            delete_if_exists(lmeXlsx); % keep clean if partial/failed
        end
    else
        fprintf(LOG, 'No successful LME outputs produced; nothing written (log only).\n');
    end
end

fclose(LOG);

%% ---------------- PCA ----------------------------------------------------
if DO_PCA_POWERS_ONLY
    run_pca_powers_only(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, phaseList, xlsxOut, dirs.Fig_PCA, ...
        SAVE_DPI, FONT_NAME, TICK_DIR, DO_PCA_CLUSTERING, PhotoperiodColourMap);
end
if DO_PCA_DELTAS_ONLY
    run_pca_deltas_only(Tpair, UR_BANDS, phaseList, xlsxOut, dirs.Fig_PCA, ...
        SAVE_DPI, FONT_NAME, TICK_DIR, DO_PCA_CLUSTERING, PhotoperiodColourMap);
end

fprintf('\nSummary-first stage complete. Outputs already written.\n');

%% ========================================================================
% PHASE B: TIME-SERIES (STREAMED)
% ========================================================================
Tdom = table();
Tcov = table();
TstabRP = table();
TstabRPOW = table();

if NEED_TS
    fprintf('\nTime-series stage starting (streamed). This is the slow bit.\n');

    stems = unique(BCS.FileStemPkg);
    stems = stems(~ismissing(stems));

    for i = 1:numel(stems)
        stem = char(stems(i));

        if ~isKey(TSMap, stem)
            fprintf('  [%d/%d] TS not mapped for %s (skipping)\n', i, numel(stems), stem);
            continue;
        end

        tsPath = TSMap(stem);
        if ~isfile(tsPath)
            fprintf('  [%d/%d] TS missing for %s: %s\n', i, numel(stems), stem, tsPath);
            continue;
        end

        fprintf('  [%d/%d] Loading TS: %s\n', i, numel(stems), tsPath);

        try
            S = load(tsPath, 'pkgTS');
            if ~isfield(S,'pkgTS') || isempty(S.pkgTS), continue; end
            pkgTS = S.pkgTS;

            BP = table(); RP = table(); RPOW = table();
            if isfield(pkgTS,'tables') && isfield(pkgTS.tables,'BandPower_Long'),    BP   = pkgTS.tables.BandPower_Long; end
            if isfield(pkgTS,'tables') && isfield(pkgTS.tables,'RidgePeriod_Long'), RP   = pkgTS.tables.RidgePeriod_Long; end
            if isfield(pkgTS,'tables') && isfield(pkgTS.tables,'RidgePower_Long'),  RPOW = pkgTS.tables.RidgePower_Long; end

            BP   = harmonise_long_table(BP);
            RP   = harmonise_long_table(RP);
            RPOW = harmonise_long_table(RPOW);

            % v7: use Raw CR and validated Raw UR only for TS analyses.
            BP   = filter_long_for_validated_raw_analysis(BP,   CarryForward_Periods, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, USE_HSUB_VALIDATION_MAP, VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK);
            RP   = filter_long_for_validated_raw_analysis(RP,   CarryForward_Periods, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, USE_HSUB_VALIDATION_MAP, VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK);
            RPOW = filter_long_for_validated_raw_analysis(RPOW, CarryForward_Periods, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, USE_HSUB_VALIDATION_MAP, VALIDATION_APPLY_ALL_PHASE_TO_LIGHT_DARK);

            % v8: stable-day filter for streamed time-series outputs. This does
            % not alter summary-level tables already produced by Script 2.
            if EXCLUDE_INITIAL_DAYS_FROM_TS_ANALYSES
                BP   = filter_time_days_minimum(BP,   MIN_DAY_FOR_STABLE_TS);
                RP   = filter_time_days_minimum(RP,   MIN_DAY_FOR_STABLE_TS);
                RPOW = filter_time_days_minimum(RPOW, MIN_DAY_FOR_STABLE_TS);
            end

            if DO_COVERAGE
                Tc = compute_true_coverage(BP, RP, RPOW, TS_PHASES);
                if ~isempty(Tc), Tcov = [Tcov; Tc]; end %#ok<AGROW>
            end

            if DO_DOMINANCE_BP && ~isempty(BP)
                Td = compute_dominance_occupancy_BP_truecov(BP, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, TS_PHASES);
                if ~isempty(Td), Tdom = [Tdom; Td]; end %#ok<AGROW>
            end

            if DO_STABILITY_RP && ~isempty(RP)
                Ts1 = compute_stability_from_long_truecov(RP, TS_PHASES, 'RidgePeriod');
                if ~isempty(Ts1), TstabRP = [TstabRP; Ts1]; end %#ok<AGROW>
            end
            if DO_STABILITY_RPOW && ~isempty(RPOW)
                Ts2 = compute_stability_from_long_truecov(RPOW, TS_PHASES, 'RidgePower');
                if ~isempty(Ts2), TstabRPOW = [TstabRPOW; Ts2]; end %#ok<AGROW>
            end

            clear pkgTS BP RP RPOW;

        catch ME
            fprintf('    TS failed for %s: %s\n', stem, ME.message);
            append_text(runLog, sprintf('TS failed for %s: %s\n', stem, ME.message));
        end
    end

    if ~isempty(Tcov),      writetable(Tcov,      xlsxOut, 'Sheet', 'Coverage_TS'); end
    if ~isempty(Tdom),      writetable(Tdom,      xlsxOut, 'Sheet', 'Dominance_BP'); end
    if ~isempty(TstabRP),   writetable(TstabRP,   xlsxOut, 'Sheet', 'Stability_RP_TS'); end
    if ~isempty(TstabRPOW), writetable(TstabRPOW, xlsxOut, 'Sheet', 'Stability_RPOW_TS'); end

    if ~isempty(Tcov)
        plot_coverage_fraction_bars(Tcov, dirs.Fig_Cov, SAVE_DPI, FONT_NAME, TICK_DIR, PhotoperiodColourMap);
    end
    if ~isempty(Tdom)
        plot_dominance_bars(Tdom, dirs.Fig_Dom, SAVE_DPI, FONT_NAME, TICK_DIR, UR_BANDS, PhotoperiodColourMap);
    end
    if ~isempty(TstabRP)
        plot_stability_bars(TstabRP, dirs.Fig_StabRP, SAVE_DPI, FONT_NAME, TICK_DIR, 'RidgePeriod', PhotoperiodColourMap);
    end
    if ~isempty(TstabRPOW)
        plot_stability_bars(TstabRPOW, dirs.Fig_StabRPOW, SAVE_DPI, FONT_NAME, TICK_DIR, 'RidgePower', PhotoperiodColourMap);
    end
end

fprintf('\nAcrossPhotoperiod analysis complete.\nOutputs:\n  %s\nWorkbook:\n  %s\n', outDir, xlsxOut);
append_text(runLog, sprintf('Run completed: %s\n', datestr(now,31)));

%% ========================================================================
% Local functions
% ========================================================================

function ensure_dir(p)
    if ~exist(p,'dir'), mkdir(p); end
end

function delete_if_exists(p)
    if exist(p,'file'), delete(p); end
end

function write_text(pathTxt, str)
    try
        fid = fopen(pathTxt, 'w');
        if fid < 0, return; end
        fprintf(fid, '%s\n', str);
        fclose(fid);
    catch
    end
end

function append_text(pathTxt, str)
    try
        fid = fopen(pathTxt, 'a');
        if fid < 0, return; end
        fprintf(fid, '%s', str);
        fclose(fid);
    catch
    end
end

function s = safe_field(nameIn)
    s = matlab.lang.makeValidName(char(string(nameIn)));
end

function sh = safe_sheet_name(nameIn)
    sh = char(string(nameIn));
    sh = regexprep(sh, '[^\w]', '_');
    if numel(sh) > 31, sh = sh(1:31); end
end


function [T, SexInfo] = add_sex_column_if_available(T, includeSex, allowInferFromSignalID)
    SexInfo = table();
    SexInfo.UseSexInLME = false;
    SexInfo.SexSource = "NotUsed";
    SexInfo.N_Male = 0;
    SexInfo.N_Female = 0;
    SexInfo.N_Unknown = height(T);

    if ~includeSex || isempty(T) || ~istable(T)
        return;
    end

    sexVar = "";
    candidates = ["Sex", "sex", "SEX", "BiologicalSex", "MouseSex"];
    for c = 1:numel(candidates)
        if ismember(candidates(c), string(T.Properties.VariableNames))
            sexVar = candidates(c);
            break;
        end
    end

    if strlength(sexVar) > 0
        rawSex = string(T.(char(sexVar)));
        SexInfo.SexSource = "ExistingColumn:" + sexVar;
    elseif allowInferFromSignalID && ismember('SignalID', T.Properties.VariableNames)
        rawSex = infer_sex_from_signalid(string(T.SignalID));
        SexInfo.SexSource = "InferredFromSignalID";
    else
        return;
    end

    sexStd = standardise_sex_labels(rawSex);
    T.Sex = categorical(sexStd);
    SexInfo.N_Male = sum(sexStd == "Male");
    SexInfo.N_Female = sum(sexStd == "Female");
    SexInfo.N_Unknown = sum(sexStd == "Unknown" | ismissing(sexStd));

    if SexInfo.N_Male > 0 && SexInfo.N_Female > 0
        SexInfo.UseSexInLME = true;
    else
        SexInfo.UseSexInLME = false;
    end
end

function tf = has_usable_sex(T)
    tf = false;
    if ~istable(T) || ~ismember('Sex', T.Properties.VariableNames)
        return;
    end
    sx = string(T.Sex);
    sx = sx(~ismissing(sx) & sx ~= "Unknown" & strlength(sx) > 0);
    tf = numel(unique(sx)) >= 2;
end

function fml = build_lme_formula(responseName, hasPhase, hasSex, sexMode)
    terms = "Photoperiod_h";
    if hasPhase
        terms = terms + " + Phase";
    end
    if hasSex
        if string(sexMode) == "PhotoperiodBySex"
            terms = terms + " + Photoperiod_h:Sex + Sex";
        else
            terms = terms + " + Sex";
        end
    end
    fml = char(string(responseName) + " ~ " + terms + " + (1|SignalID) + (1|File)");
end

function sx = infer_sex_from_signalid(signalID)
    s = lower(string(signalID));
    sx = repmat("Unknown", numel(s), 1);
    % Conservative patterns: explicit words or clear boundary-delimited M/F labels.
    isFemale = contains(s, "female") | contains(s, "_f_") | contains(s, "-f-") | ...
               startsWith(s, "f_") | startsWith(s, "f-") | endsWith(s, "_f") | endsWith(s, "-f") | ...
               ~cellfun('isempty', regexp(cellstr(s), '(^|[^a-z])f[0-9]+($|[^a-z0-9])', 'once'));
    isMale = contains(s, "male") | contains(s, "_m_") | contains(s, "-m-") | ...
             startsWith(s, "m_") | startsWith(s, "m-") | endsWith(s, "_m") | endsWith(s, "-m") | ...
             ~cellfun('isempty', regexp(cellstr(s), '(^|[^a-z])m[0-9]+($|[^a-z0-9])', 'once'));
    % Avoid classifying 'female' as male because it contains 'male'.
    isMale = isMale & ~contains(s, "female");
    sx(isMale) = "Male";
    sx(isFemale) = "Female";
end

function sx = standardise_sex_labels(rawSex)
    s = lower(strtrim(string(rawSex)));
    sx = repmat("Unknown", numel(s), 1);
    sx(s == "m" | s == "male" | s == "man" | s == "1") = "Male";
    sx(s == "f" | s == "female" | s == "woman" | s == "2") = "Female";
    sx(contains(s,"female")) = "Female";
    sx(contains(s,"male") & ~contains(s,"female")) = "Male";
end

function Tinf = annotate_lme_inference_table(T, fdrFamily, metricClass, bandName, responseName, formulaStr, tableType)
    % Converts MATLAB fitlme coefficient/ANOVA tables into a harmonised
    % inference table suitable for transparent BH/FDR correction. Intercepts
    % are retained but excluded from FDR because they are not biological
    % hypothesis tests.
    if isempty(T) || ~istable(T)
        Tinf = table();
        return;
    end

    n = height(T);
    term = get_lme_terms(T);
    pRaw = get_numeric_var_or_nan(T, {'pValue','PValue','pVal','Prob_F','ProbF'});

    isIntercept = contains(lower(term), 'intercept');
    include = isfinite(pRaw) & ~isIntercept;
    if string(tableType) == "ANOVA"
        include = isfinite(pRaw);
    end

    Tinf = table();
    Tinf.FDRFamily      = repmat(string(fdrFamily), n, 1);
    Tinf.MetricClass    = repmat(string(metricClass), n, 1);
    Tinf.BandName       = repmat(string(bandName), n, 1);
    Tinf.Response       = repmat(string(responseName), n, 1);
    Tinf.TableType      = repmat(string(tableType), n, 1);
    Tinf.UR12_18SensitivityFlag = repmat(contains(string(bandName), 'UR_12_18'), n, 1);
    Tinf.Term           = term(:);
    Tinf.Formula        = repmat(string(formulaStr), n, 1);
    Tinf.p_raw          = pRaw(:);
    Tinf.IncludeInFDR   = include(:);

    % Common model statistics, if present. Missing quantities stay NaN.
    Tinf.Estimate = get_numeric_var_or_nan(T, {'Estimate'});
    Tinf.SE       = get_numeric_var_or_nan(T, {'SE'});
    Tinf.tStat    = get_numeric_var_or_nan(T, {'tStat'});
    Tinf.FStat    = get_numeric_var_or_nan(T, {'FStat'});
    Tinf.DF       = get_numeric_var_or_nan(T, {'DF'});
    Tinf.DF1      = get_numeric_var_or_nan(T, {'DF1'});
    Tinf.DF2      = get_numeric_var_or_nan(T, {'DF2'});
end

function term = get_lme_terms(T)
    n = height(T);
    term = strings(n,1);
    if ismember('Name', T.Properties.VariableNames)
        term = string(T.Name);
        return;
    end
    if ~isempty(T.Properties.RowNames) && numel(T.Properties.RowNames) == n
        term = string(T.Properties.RowNames(:));
        return;
    end
    if ismember('Term', T.Properties.VariableNames)
        term = string(T.Term);
        return;
    end
    for i = 1:n
        term(i) = "Row_" + string(i);
    end
end

function x = get_numeric_var_or_nan(T, nameCandidates)
    n = height(T);
    x = nan(n,1);
    for k = 1:numel(nameCandidates)
        nm = char(nameCandidates{k});
        if ismember(nm, T.Properties.VariableNames)
            try
                v = T.(nm);
                if iscell(v), v = string(v); end
                x = double(v);
                x = x(:);
                if numel(x) ~= n, x = nan(n,1); end
                return;
            catch
                x = nan(n,1);
                return;
            end
        end
    end
end

function Tout = apply_bh_fdr_by_family(Tin, alphaVal)
    Tout = Tin;
    n = height(Tout);
    Tout.p_BH = nan(n,1);
    Tout.Significant_BH = false(n,1);
    Tout.FDR_Rank = nan(n,1);
    Tout.FDR_m = nan(n,1);
    Tout.FDR_Alpha = repmat(alphaVal, n, 1);

    if isempty(Tout) || ~ismember('FDRFamily', Tout.Properties.VariableNames) || ~ismember('p_raw', Tout.Properties.VariableNames)
        return;
    end
    if ~ismember('IncludeInFDR', Tout.Properties.VariableNames)
        Tout.IncludeInFDR = isfinite(Tout.p_raw);
    end

    fams = unique(string(Tout.FDRFamily));
    fams = fams(~ismissing(fams));
    for f = 1:numel(fams)
        idxFam = string(Tout.FDRFamily) == fams(f) & Tout.IncludeInFDR & isfinite(Tout.p_raw);
        ii = find(idxFam);
        if isempty(ii), continue; end
        p = Tout.p_raw(ii);
        [q, rankOriginal] = bh_adjust(p);
        Tout.p_BH(ii) = q;
        Tout.Significant_BH(ii) = q <= alphaVal;
        Tout.FDR_Rank(ii) = rankOriginal;
        Tout.FDR_m(ii) = numel(ii);
    end
end

function [q, rankOriginal] = bh_adjust(p)
    p = p(:);
    m = numel(p);
    q = nan(m,1);
    rankOriginal = nan(m,1);
    if m == 0, return; end
    [ps, ord] = sort(p, 'ascend');
    ranks = (1:m)';
    qs = ps .* m ./ ranks;
    % Enforce monotonicity from largest to smallest p-value.
    for i = m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end
    qs = min(qs, 1);
    q(ord) = qs;
    rankOriginal(ord) = ranks;
end

function write_lme_fdr_subset(T, xlsxPath, metricClass, tableType, sheetName)
    try
        if isempty(T) || ~istable(T), return; end
        idx = string(T.MetricClass) == string(metricClass) & string(T.TableType) == string(tableType);
        S = T(idx, :);
        if ~isempty(S)
            writetable(S, xlsxPath, 'Sheet', safe_sheet_name(sheetName));
        end
    catch
    end
end

function A = append_table_rows(A, B)
    if isempty(B) || ~istable(B)
        return;
    end
    if isempty(A) || ~istable(A) || width(A) == 0
        A = B;
    else
        A = [A; B]; %#ok<AGROW>
    end
end

function out = getfield_or(S, fld, defaultVal)
    out = defaultVal;
    try
        if isstruct(S) && isfield(S, fld)
            out = S.(fld);
        end
    catch
        out = defaultVal;
    end
end

function stem = extract_stem_from_summary(summaryMatPath)
    stem = '';
    try
        [~, nm, ~] = fileparts(summaryMatPath);
        nm = string(nm);
        pref = "WP_Summary__";
        if startsWith(nm, pref)
            stem = char(extractAfter(nm, strlength(pref)));
        end
    catch
        stem = '';
    end
end

function T = harmonise_long_table(T)
    if isempty(T), return; end
    if ismember('File', T.Properties.VariableNames), T.File = string(T.File); end
    if ismember('SignalID', T.Properties.VariableNames), T.SignalID = string(T.SignalID); end
    if ismember('Source', T.Properties.VariableNames), T.Source = string(T.Source); end
    if ismember('BandName', T.Properties.VariableNames), T.BandName = string(T.BandName); end
    if ismember('Phase', T.Properties.VariableNames), T.Phase = string(T.Phase); end
    if ismember('LightStateValue', T.Properties.VariableNames), T.LightStateValue = string(T.LightStateValue); end
end

function apply_axes_style(ax, FONT_NAME, TICK_DIR)
    if isempty(ax) || ~isvalid(ax), return; end
    ax.FontName = FONT_NAME;
    ax.TickDir  = TICK_DIR;
    ax.Box      = 'off';      % removes top/right axes
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.LineWidth = 1.1;
    ax.TickLength = [0.012 0.012];
end

function set_bold_labels(ttl, xl, yl)
    try
        if ~isempty(ttl), set(ttl,'FontWeight','bold'); end
        if ~isempty(xl),  set(xl,'FontWeight','bold');  end
        if ~isempty(yl),  set(yl,'FontWeight','bold');  end
    catch
    end
end

function print_jpeg600(fig, outFn, SAVE_DPI)
    try
        print(fig, outFn, '-djpeg', sprintf('-r%d', SAVE_DPI));
    catch
    end
end

function s = sanitise_filename(strIn)
    s = regexprep(string(strIn), '[^\w\-]', '_');
    s = char(s);
    if isempty(s), s = 'x'; end
    if numel(s) > 160, s = s(1:160); end
end

%% ---------------- COLOUR GOVERNANCE -------------------------------------

function base = tol_bright_base()
    % Paul Tol "Bright" qualitative palette (7 colours)
    base = [
        68 119 170;
        102 204 238;
        34 136 51;
        204 187 68;
        238 102 119;
        170 51 119;
        187 187 187
    ] / 255;
end

function cmap = tol_bright(n)
    base = tol_bright_base();
    if n <= size(base,1)
        cmap = base(1:n,:);
    else
        cmap = lines(n);
        cmap(1:size(base,1),:) = base;
    end
end

function [bandMap, ppMap] = try_get_colour_maps_from_pkgS(pkgS)
    % Robust attempt to recover colour maps from Script 1 handoff.
    % Accept either:
    %   - pkgS.colors.BandColourMap / pkgS.colors.PhotoperiodColourMap
    %   - pkgS.meta.colors.* variants
    %   - pkgS.ColorMaps.* variants
    bandMap = [];
    ppMap = [];

    cands = {
        {'colors'}
        {'meta','colors'}
        {'ColorMaps'}
        {'meta','ColorMaps'}
    };

    for i = 1:numel(cands)
        node = get_nested_field(pkgS, cands{i});
        if isempty(node), continue; end

        try
            if isstruct(node)
                if isfield(node,'BandColourMap'), bandMap = node.BandColourMap; end
                if isfield(node,'BandColorMap') && isempty(bandMap), bandMap = node.BandColorMap; end

                if isfield(node,'PhotoperiodColourMap'), ppMap = node.PhotoperiodColourMap; end
                if isfield(node,'PhotoperiodColorMap') && isempty(ppMap), ppMap = node.PhotoperiodColorMap; end

                % Sometimes stored as table
                if isfield(node,'BandColoursTable') && isempty(bandMap), bandMap = node.BandColoursTable; end
                if isfield(node,'PhotoperiodColoursTable') && isempty(ppMap), ppMap = node.PhotoperiodColoursTable; end
            end
        catch
        end

        if ~isempty(bandMap) || ~isempty(ppMap)
            return;
        end
    end
end

function node = get_nested_field(S, pathCell)
    node = S;
    for k = 1:numel(pathCell)
        key = pathCell{k};
        if isstruct(node) && isfield(node, key)
            node = node.(key);
        else
            node = [];
            return;
        end
    end
end

function BandColourMap = build_band_colour_map(BANDS_ALL, handoffObj)
    % Returns containers.Map band -> 1x3 colour row
    BandColourMap = containers.Map('KeyType','char','ValueType','any');

    % 1) Try handoff
    if ~isempty(handoffObj)
        try
            BandColourMap = parse_colour_map_any(handoffObj, "BandName", "Band", "Component");
            if validate_map_has_keys(BandColourMap, cellstr(BANDS_ALL))
                return;
            end
        catch
        end
    end

    % 2) Fallback: assign Tol Bright uniquely to each band in fixed order
    cols = tol_bright(numel(BANDS_ALL));
    for i = 1:numel(BANDS_ALL)
        BandColourMap(char(BANDS_ALL(i))) = cols(i,:);
    end
end

function PhotoperiodColourMap = build_photoperiod_colour_map(ppValsAll, handoffObj)
    % Returns containers.Map photoperiod numeric-as-string -> 1x3 colour row
    PhotoperiodColourMap = containers.Map('KeyType','char','ValueType','any');

    % 1) Try handoff
    if ~isempty(handoffObj)
        try
            PhotoperiodColourMap = parse_colour_map_any(handoffObj, "Photoperiod_h", "Photoperiod", "PP");
            keysNeed = cellstr(string(ppValsAll));
            if validate_map_has_keys(PhotoperiodColourMap, keysNeed)
                return;
            end
        catch
        end
    end

    % 2) Fallback: Tol Bright by sorted photoperiods
    cols = tol_bright(numel(ppValsAll));
    for i = 1:numel(ppValsAll)
        PhotoperiodColourMap(char(string(ppValsAll(i)))) = cols(i,:);
    end
end

function ok = validate_map_has_keys(M, keysNeed)
    ok = true;
    for i = 1:numel(keysNeed)
        if ~isKey(M, char(keysNeed{i}))
            ok = false;
            return;
        end
    end
end

function M = parse_colour_map_any(obj, varargin)
    % Accept:
    %   - containers.Map with char keys -> RGB rows
    %   - struct with fields named keys and 1x3 values
    %   - table with key column + RGB columns (R,G,B or RGB)
    M = containers.Map('KeyType','char','ValueType','any');

    if isa(obj,'containers.Map')
        k = obj.keys;
        for i = 1:numel(k)
            key = char(string(k{i}));
            v = obj(k{i});
            v = normalise_rgb(v);
            if ~isempty(v), M(key) = v; end
        end
        return;
    end

    if isstruct(obj)
        fn = fieldnames(obj);
        for i = 1:numel(fn)
            key = char(string(fn{i}));
            v = obj.(fn{i});
            v = normalise_rgb(v);
            if ~isempty(v), M(key) = v; end
        end
        return;
    end

    if istable(obj)
        % Find key column
        keyCol = '';
        for i = 1:numel(varargin)
            nm = char(string(varargin{i}));
            if ismember(nm, obj.Properties.VariableNames)
                keyCol = nm;
                break;
            end
        end
        if isempty(keyCol)
            keyCol = obj.Properties.VariableNames{1};
        end

        % RGB columns
        if all(ismember({'R','G','B'}, obj.Properties.VariableNames))
            R = obj.R; G = obj.G; B = obj.B;
            for i = 1:height(obj)
                key = char(string(obj.(keyCol)(i)));
                v = [double(R(i)) double(G(i)) double(B(i))];
                v = normalise_rgb(v);
                if ~isempty(v), M(key) = v; end
            end
            return;
        end

        if ismember('RGB', obj.Properties.VariableNames)
            for i = 1:height(obj)
                key = char(string(obj.(keyCol)(i)));
                v = obj.RGB(i,:);
                v = normalise_rgb(v);
                if ~isempty(v), M(key) = v; end
            end
            return;
        end

        % Try three numeric columns after key
        try
            vi = find(strcmp(obj.Properties.VariableNames, keyCol), 1, 'first');
            numCols = obj.Properties.VariableNames;
            numCols(vi) = [];
            numCols = numCols(1:min(3,numel(numCols)));
            if numel(numCols) >= 3
                for i = 1:height(obj)
                    key = char(string(obj.(keyCol)(i)));
                    v = [double(obj.(numCols{1})(i)) double(obj.(numCols{2})(i)) double(obj.(numCols{3})(i))];
                    v = normalise_rgb(v);
                    if ~isempty(v), M(key) = v; end
                end
            end
        catch
        end
    end
end

function v = normalise_rgb(vIn)
    v = [];
    try
        vIn = double(vIn);
        if isvector(vIn) && numel(vIn) == 3
            v = vIn(:).';
            if any(v > 1.5) % assume 0..255
                v = v ./ 255;
            end
            v = max(0, min(1, v));
        end
    catch
        v = [];
    end
end

%% ---------------- SUMMARY TABLE BUILDERS --------------------------------

function [Tpair, Tsum] = build_CR_UR_pairs_from_BCS(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    CR = BCS(BCS.Source==SRC_CR & BCS.BandName==CR_BAND, :);
    CR = CR(:, {'File','SignalID','Photoperiod_h','Phase','MeanBandPower_log10','FileStemPkg'});
    CR.Properties.VariableNames{'MeanBandPower_log10'} = 'CR_Log10';

    UR = BCS(BCS.Source==SRC_UR & ismember(BCS.BandName, UR_BANDS), :);
    UR = UR(:, {'File','SignalID','Photoperiod_h','Phase','BandName','MeanBandPower_log10','FileStemPkg'});
    UR.Properties.VariableNames{'BandName'} = 'UR_Band';
    UR.Properties.VariableNames{'MeanBandPower_log10'} = 'UR_Log10';

    keys = {'File','SignalID','Photoperiod_h','Phase','FileStemPkg'};
    Tpair = innerjoin(UR, CR, 'Keys', keys);

    Tpair.Delta_log10 = Tpair.UR_Log10 - Tpair.CR_Log10;
    Tpair.Ratio = 10 .^ (Tpair.Delta_log10);

    G = findgroups(Tpair.Photoperiod_h, Tpair.Phase, Tpair.UR_Band);
    meanCR = splitapply(@(x) mean(x,'omitnan'), Tpair.CR_Log10, G);
    sdCR   = splitapply(@(x) std(x,'omitnan'),  Tpair.CR_Log10, G);
    meanUR = splitapply(@(x) mean(x,'omitnan'), Tpair.UR_Log10, G);
    sdUR   = splitapply(@(x) std(x,'omitnan'),  Tpair.UR_Log10, G);
    meanD  = splitapply(@(x) mean(x,'omitnan'), Tpair.Delta_log10, G);
    sdD    = splitapply(@(x) std(x,'omitnan'),  Tpair.Delta_log10, G);

    [pp, ph, ub] = splitapply(@(a,b,c) deal(a(1),b(1),c(1)), Tpair.Photoperiod_h, Tpair.Phase, Tpair.UR_Band, G);

    Tsum = table(pp, ph, ub, meanCR, sdCR, meanUR, sdUR, meanD, sdD, ...
        'VariableNames', {'Photoperiod_h','Phase','UR_Band', ...
                          'Mean_CR_Log10','SD_CR_Log10','Mean_UR_Log10','SD_UR_Log10', ...
                          'Mean_Delta_log10','SD_Delta_log10'});
    Tsum = sortrows(Tsum, {'Phase','UR_Band','Photoperiod_h'});
end

function Tabs = build_absolute_summary_tables(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    bandsAll = [CR_BAND UR_BANDS];
    rows = cell(0,5);
    for b = 1:numel(bandsAll)
        bn = bandsAll(b);
        if bn == CR_BAND
            D = BCS(BCS.Source==SRC_CR & BCS.BandName==bn, :);
        else
            D = BCS(BCS.Source==SRC_UR & BCS.BandName==bn, :);
        end
        if isempty(D), continue; end

        G = findgroups(D.Photoperiod_h, D.Phase);
        pp = splitapply(@(x) x(1), D.Photoperiod_h, G);
        ph = splitapply(@(x) x(1), D.Phase, G);
        m  = splitapply(@(x) mean(x,'omitnan'), D.MeanBandPower_log10, G);
        s  = splitapply(@(x) std(x,'omitnan'),  D.MeanBandPower_log10, G);

        n = numel(m);
        add = [num2cell(pp(:)) cellstr(ph(:)) cellstr(repmat(string(bn), n, 1)) num2cell(m(:)) num2cell(s(:))];
        rows = [rows; add]; %#ok<AGROW>
    end

    Tabs = cell2table(rows, 'VariableNames', {'Photoperiod_h','Phase','BandName','Mean_Log10','SD_Log10'});
    Tabs.Phase = string(Tabs.Phase);
    Tabs.BandName = string(Tabs.BandName);
    Tabs = sortrows(Tabs, {'Phase','BandName','Photoperiod_h'});
end

function TURcomp = build_UR_composition_linearfrac(BCS, UR_BANDS, SRC_UR)
    D = BCS(BCS.Source==SRC_UR & ismember(BCS.BandName, UR_BANDS), :);
    if isempty(D), TURcomp = table(); return; end

    keys = {'File','SignalID','Photoperiod_h','Phase'};
    Ukeys = unique(D(:, keys));

    rows = cell(0,4);
    for i = 1:height(Ukeys)
        f  = Ukeys.File(i);
        id = Ukeys.SignalID(i);
        pp = Ukeys.Photoperiod_h(i);
        ph = Ukeys.Phase(i);

        d0 = D(D.File==f & D.SignalID==id & D.Photoperiod_h==pp & D.Phase==ph, :);
        if isempty(d0), continue; end

        linP = NaN(numel(UR_BANDS),1);
        for b = 1:numel(UR_BANDS)
            r = d0(d0.BandName==UR_BANDS(b), :);
            if ~isempty(r)
                linP(b) = 10 .^ r.MeanBandPower_log10(1);
            end
        end
        denom = sum(linP, 'omitnan');
        if isfinite(denom) && denom > 0
            frac = linP ./ denom;
        else
            frac = NaN(numel(UR_BANDS),1);
        end

        n = numel(UR_BANDS);
        add = [num2cell(repmat(pp, n, 1)) cellstr(repmat(string(ph), n, 1)) cellstr(UR_BANDS(:)) num2cell(frac(:))];
        rows = [rows; add]; %#ok<AGROW>
    end

    T = cell2table(rows, 'VariableNames', {'Photoperiod_h','Phase','UR_Band','Frac_linear'});
    T.Phase = string(T.Phase);
    T.UR_Band = string(T.UR_Band);

    G = findgroups(T.Photoperiod_h, T.Phase, T.UR_Band);
    m = splitapply(@(x) mean(x,'omitnan'), T.Frac_linear, G);
    s = splitapply(@(x) std(x,'omitnan'),  T.Frac_linear, G);
    [pp, ph, ub] = splitapply(@(a,b,c) deal(a(1),b(1),c(1)), T.Photoperiod_h, T.Phase, T.UR_Band, G);

    TURcomp = table(pp, ph, ub, m, s, 'VariableNames', {'Photoperiod_h','Phase','UR_Band','MeanFrac','SDFrac'});
    TURcomp = sortrows(TURcomp, {'Phase','UR_Band','Photoperiod_h'});
end

function TCompAll = build_CRplusUR_composition_linearfrac(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    % Fractions across [CR + all UR bands], linear power domain.
    keys = {'File','SignalID','Photoperiod_h','Phase'};

    CR = BCS(BCS.Source==SRC_CR & BCS.BandName==CR_BAND, :);
    UR = BCS(BCS.Source==SRC_UR & ismember(BCS.BandName, UR_BANDS), :);

    if isempty(CR) || isempty(UR)
        TCompAll = table();
        return;
    end

    Ukeys = unique([CR(:,keys); UR(:,keys)]);
    rows = cell(0,4);

    for i = 1:height(Ukeys)
        f  = Ukeys.File(i);
        id = Ukeys.SignalID(i);
        pp = Ukeys.Photoperiod_h(i);
        ph = Ukeys.Phase(i);

        pCR = NaN;
        rCR = CR(CR.File==f & CR.SignalID==id & CR.Photoperiod_h==pp & CR.Phase==ph, :);
        if ~isempty(rCR)
            pCR = 10 .^ rCR.MeanBandPower_log10(1);
        end

        pUR = NaN(numel(UR_BANDS),1);
        rUR = UR(UR.File==f & UR.SignalID==id & UR.Photoperiod_h==pp & UR.Phase==ph, :);

        for b = 1:numel(UR_BANDS)
            rr = rUR(rUR.BandName==UR_BANDS(b), :);
            if ~isempty(rr)
                pUR(b) = 10 .^ rr.MeanBandPower_log10(1);
            end
        end

        denom = pCR + sum(pUR,'omitnan');
        if isfinite(denom) && denom > 0
            fracCR = pCR ./ denom;
            fracUR = pUR ./ denom;
        else
            fracCR = NaN;
            fracUR = NaN(numel(UR_BANDS),1);
        end

        % write CR row
        rows(end+1,:) = {pp, ph, "CR", fracCR}; %#ok<AGROW>
        % write UR rows
        for b = 1:numel(UR_BANDS)
            rows(end+1,:) = {pp, ph, UR_BANDS(b), fracUR(b)}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', {'Photoperiod_h','Phase','Component','Frac_linear'});
    T.Phase = string(T.Phase);
    T.Component = string(T.Component);

    G = findgroups(T.Photoperiod_h, T.Phase, T.Component);
    m = splitapply(@(x) mean(x,'omitnan'), T.Frac_linear, G);
    s = splitapply(@(x) std(x,'omitnan'),  T.Frac_linear, G);
    [pp, ph, comp] = splitapply(@(a,b,c) deal(a(1),b(1),c(1)), T.Photoperiod_h, T.Phase, T.Component, G);

    TCompAll = table(pp, ph, comp, m, s, 'VariableNames', {'Photoperiod_h','Phase','Component','MeanFrac','SDFrac'});
    TCompAll = sortrows(TCompAll, {'Phase','Component','Photoperiod_h'});
end

function TBP = build_bandpower_table_from_BCS(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    rows = cell(0,7);
    hdr = {'File','SignalID','Photoperiod_h','Phase','BandName','Source','MeanBandPower_log10'};

    CR = BCS(BCS.Source==SRC_CR & BCS.BandName==CR_BAND, :);
    if ~isempty(CR)
        add = [cellstr(CR.File) cellstr(CR.SignalID) num2cell(CR.Photoperiod_h) cellstr(CR.Phase) cellstr(CR.BandName) cellstr(CR.Source) num2cell(CR.MeanBandPower_log10)];
        rows = [rows; add]; %#ok<AGROW>
    end

    UR = BCS(BCS.Source==SRC_UR & ismember(BCS.BandName, UR_BANDS), :);
    if ~isempty(UR)
        add = [cellstr(UR.File) cellstr(UR.SignalID) num2cell(UR.Photoperiod_h) cellstr(UR.Phase) cellstr(UR.BandName) cellstr(UR.Source) num2cell(UR.MeanBandPower_log10)];
        rows = [rows; add]; %#ok<AGROW>
    end

    TBP = cell2table(rows, 'VariableNames', hdr);
    TBP.File = string(TBP.File);
    TBP.SignalID = string(TBP.SignalID);
    TBP.Phase = string(TBP.Phase);
    TBP.BandName = string(TBP.BandName);
    TBP.Source = string(TBP.Source);
end

%% ---------------- SUMMARY PLOT HELPERS ----------------------------------

function plot_absolute_band_across_pp(Tband, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, titleStr, PPColourMap)
    if isempty(Tband), return; end
    phases = unique(Tband.Phase);

    for ph = phases(:).'
        D = Tband(Tband.Phase==ph, :);
        if isempty(D), continue; end

        fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
        ax = gca;
        apply_axes_style(ax, FONT_NAME, TICK_DIR);

        pp = double(D.Photoperiod_h);
        ppU = unique(pp);
        ppU = sort(ppU(isfinite(ppU)));

        m = NaN(numel(ppU),1);
        s = NaN(numel(ppU),1);
        for i = 1:numel(ppU)
            r = D(pp==ppU(i), :);
            if ~isempty(r)
                m(i) = r.Mean_Log10(1);
                s(i) = r.SD_Log10(1);
            end
        end

        X = categorical(string(ppU));
        bb = bar(X, m, 'FaceColor','flat');
        for i = 1:numel(ppU)
            key = char(string(ppU(i)));
            if isKey(PPColourMap, key)
                bb.CData(i,:) = PPColourMap(key);
            end
        end

        hold on;
        errorbar(bb.XEndPoints, m, s, 'k', 'LineStyle','none', 'LineWidth', 1.2);
        hold off;

        yl = ylabel('Mean band power (log_{10})');
        xl = xlabel('Photoperiod (h)');
        ttl = title(sprintf('%s | Phase: %s', titleStr, char(ph)), 'Interpreter','none');
        set_bold_labels(ttl, xl, yl);

        outFn = fullfile(outDir, sprintf('Abs_%s_Phase_%s.jpg', sanitise_filename(char(D.BandName(1))), sanitise_filename(char(ph))));
        print_jpeg600(fig, outFn, SAVE_DPI);
        close(fig);
    end
end

function S = summarise_CR_UR(D)
    G = findgroups(D.Photoperiod_h);
    pp = splitapply(@(x) x(1), D.Photoperiod_h, G);
    mCR = splitapply(@(x) mean(x,'omitnan'), D.CR_Log10, G);
    sCR = splitapply(@(x) std(x,'omitnan'),  D.CR_Log10, G);
    mUR = splitapply(@(x) mean(x,'omitnan'), D.UR_Log10, G);
    sUR = splitapply(@(x) std(x,'omitnan'),  D.UR_Log10, G);
    S = table(pp, mCR, sCR, mUR, sUR, 'VariableNames', {'Photoperiod_h','Mean_CR_Log10','SD_CR_Log10','Mean_UR_Log10','SD_UR_Log10'});
    S = sortrows(S,'Photoperiod_h');
end

function Sd = summarise_delta(D)
    G = findgroups(D.Photoperiod_h);
    pp = splitapply(@(x) x(1), D.Photoperiod_h, G);
    m = splitapply(@(x) mean(x,'omitnan'), D.Delta_log10, G);
    s = splitapply(@(x) std(x,'omitnan'),  D.Delta_log10, G);
    Sd = table(pp, m, s, 'VariableNames', {'Photoperiod_h','MeanDelta','SDDelta'});
    Sd = sortrows(Sd,'Photoperiod_h');
end

function Sr = summarise_ratio(D, robustFlag)
    G = findgroups(D.Photoperiod_h);
    pp = splitapply(@(x) x(1), D.Photoperiod_h, G);

    if robustFlag
        med = splitapply(@(x) median(x,'omitnan'), D.Ratio, G);
        q1  = splitapply(@(x) prctile(x,25), D.Ratio, G);
        q3  = splitapply(@(x) prctile(x,75), D.Ratio, G);
        Sr = table(pp, med, q1, q3, 'VariableNames', {'Photoperiod_h','MedianRatio','Q1','Q3'});
    else
        m = splitapply(@(x) mean(x,'omitnan'), D.Ratio, G);
        s = splitapply(@(x) std(x,'omitnan'),  D.Ratio, G);
        Sr = table(pp, m, s, 'VariableNames', {'Photoperiod_h','MeanRatio','SDRatio'});
    end
    Sr = sortrows(Sr,'Photoperiod_h');
end

function plot_CR_vs_UR(S, CR_BAND, ub, phaseTag, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, BandColourMap)
    fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
    ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

    X = categorical(string(S.Photoperiod_h));
    Y = [S.Mean_CR_Log10, S.Mean_UR_Log10];
    bb = bar(X, Y, 'FaceColor','flat');

    % Apply band colours (consistent across bandwise outputs)
    cCR = [0 0 0];
    cUR = [0 0 0];
    if isKey(BandColourMap, char(CR_BAND)), cCR = BandColourMap(char(CR_BAND)); end
    if isKey(BandColourMap, char(ub)),      cUR = BandColourMap(char(ub)); end
    bb(1).CData = repmat(cCR, numel(S.Photoperiod_h), 1);
    bb(2).CData = repmat(cUR, numel(S.Photoperiod_h), 1);

    hold on;
    errorbar(bb(1).XEndPoints, S.Mean_CR_Log10, S.SD_CR_Log10, 'k', 'LineStyle','none', 'LineWidth', 1.2);
    errorbar(bb(2).XEndPoints, S.Mean_UR_Log10, S.SD_UR_Log10, 'k', 'LineStyle','none', 'LineWidth', 1.2);
    hold off;

    yl = ylabel('Mean band power (log_{10})');
    xl = xlabel('Photoperiod (h)');
    ttl = title(sprintf('CR vs %s across photoperiod | Phase: %s', char(ub), phaseTag), 'Interpreter','none');
    set_bold_labels(ttl, xl, yl);

    lg = legend({sprintf('%s (Raw)', char(CR_BAND)), sprintf('%s (validated Raw UR)', char(ub))}, 'Location','eastoutside','Interpreter','none');
    set(lg,'Box','off','FontName',FONT_NAME);

    outFn = fullfile(outDir, sprintf('CR_vs_%s_Phase_%s.jpg', sanitise_filename(char(ub)), sanitise_filename(phaseTag)));
    print_jpeg600(fig, outFn, SAVE_DPI);
    close(fig);
end

function plot_delta(Sd, ub, phaseTag, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, PPColourMap)
    fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
    ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

    pp = double(Sd.Photoperiod_h);
    [pp, ix] = sort(pp);
    m = Sd.MeanDelta(ix);
    s = Sd.SDDelta(ix);

    X = categorical(string(pp));
    bb = bar(X, m, 'FaceColor','flat');
    for i = 1:numel(pp)
        key = char(string(pp(i)));
        if isKey(PPColourMap, key)
            bb.CData(i,:) = PPColourMap(key);
        end
    end

    hold on;
    errorbar(bb.XEndPoints, m, s, 'k', 'LineStyle','none', 'LineWidth', 1.2);
    yline(0,'k-','LineWidth',1.0);
    hold off;

    yl = ylabel('log_{10}(UR/CR)  (UR_{band} - CR)');
    xl = xlabel('Photoperiod (h)');
    ttl = title(sprintf('Delta log10 (UR-CR): %s | Phase: %s', char(ub), phaseTag), 'Interpreter','none');
    set_bold_labels(ttl, xl, yl);

    outFn = fullfile(outDir, sprintf('DeltaLog10_URminusCR_%s_Phase_%s.jpg', sanitise_filename(char(ub)), sanitise_filename(phaseTag)));
    print_jpeg600(fig, outFn, SAVE_DPI);
    close(fig);
end

function plot_ratio(Sr, ub, phaseTag, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, robustFlag, PPColourMap)
    fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
    ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

    pp = double(Sr.Photoperiod_h);
    [pp, ix] = sort(pp);

    X = categorical(string(pp));

    if robustFlag
        med = Sr.MedianRatio(ix);
        q1  = Sr.Q1(ix);
        q3  = Sr.Q3(ix);

        bb = bar(X, med, 'FaceColor','flat');
        for i = 1:numel(pp)
            key = char(string(pp(i)));
            if isKey(PPColourMap, key)
                bb.CData(i,:) = PPColourMap(key);
            end
        end

        hold on;
        lo = med - q1;
        hi = q3 - med;
        errorbar(bb.XEndPoints, med, lo, hi, 'k', 'LineStyle','none', 'LineWidth', 1.2);
        yline(1,'k-','LineWidth',1.0);
        hold off;

        yl = ylabel('UR/CR ratio (median, IQR)');
    else
        mn = Sr.MeanRatio(ix);
        sd = Sr.SDRatio(ix);

        bb = bar(X, mn, 'FaceColor','flat');
        for i = 1:numel(pp)
            key = char(string(pp(i)));
            if isKey(PPColourMap, key)
                bb.CData(i,:) = PPColourMap(key);
            end
        end

        hold on;
        errorbar(bb.XEndPoints, mn, sd, 'k', 'LineStyle','none', 'LineWidth', 1.2);
        yline(1,'k-','LineWidth',1.0);
        hold off;

        yl = ylabel('UR/CR ratio (mean ± SD)');
    end

    xl = xlabel('Photoperiod (h)');
    ttl = title(sprintf('UR/CR ratio: %s | Phase: %s', char(ub), phaseTag), 'Interpreter','none');
    set_bold_labels(ttl, xl, yl);

    outFn = fullfile(outDir, sprintf('Ratio_URoverCR_%s_Phase_%s.jpg', sanitise_filename(char(ub)), sanitise_filename(phaseTag)));
    print_jpeg600(fig, outFn, SAVE_DPI);
    close(fig);
end

function plot_UR_composition(TURcomp, UR_BANDS, phaseList, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, BandColourMap)
    if isempty(TURcomp), return; end
    for ph = phaseList(:).'
        phaseTag = char(ph);
        D = TURcomp(TURcomp.Phase==string(phaseTag), :);
        if isempty(D), continue; end

        ppVals = unique(D.Photoperiod_h);
        ppVals = sort(ppVals(~isnan(ppVals)));
        MeanFrac = NaN(numel(ppVals), numel(UR_BANDS));
        SDFrac   = NaN(size(MeanFrac));

        for i = 1:numel(ppVals)
            dpp = D(D.Photoperiod_h==ppVals(i), :);
            for b = 1:numel(UR_BANDS)
                ub = UR_BANDS(b);
                r = dpp(dpp.UR_Band==ub, :);
                if ~isempty(r)
                    MeanFrac(i,b) = r.MeanFrac(1);
                    SDFrac(i,b)   = r.SDFrac(1);
                end
            end
        end

        fig = figure('Visible','off'); set(fig,'Position',[100 100 1400 560]);
        ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

        X = categorical(string(ppVals));
        bb = bar(X, MeanFrac, 'FaceColor','flat');

        for b = 1:numel(UR_BANDS)
            ub = UR_BANDS(b);
            c = [0 0 0];
            if isKey(BandColourMap, char(ub)), c = BandColourMap(char(ub)); end
            bb(b).CData = repmat(c, numel(ppVals), 1);
        end

        hold on;
        for b = 1:numel(UR_BANDS)
            errorbar(bb(b).XEndPoints, MeanFrac(:,b), SDFrac(:,b), 'k', 'LineStyle','none', 'LineWidth', 1.0);
        end
        hold off;

        ylim([0 1]);
        yl = ylabel('Fraction within UR (linear)');
        xl = xlabel('Photoperiod (h)');
        ttl = title(sprintf('UR composition (linear fractions) | Phase: %s', phaseTag), 'Interpreter','none');
        set_bold_labels(ttl, xl, yl);

        lg = legend(cellstr(UR_BANDS), 'Location','eastoutside', 'Interpreter','none');
        set(lg,'Box','off','FontName',FONT_NAME);

        outFn = fullfile(outDir, sprintf('UR_FracCompositionLinear_Phase_%s.jpg', sanitise_filename(phaseTag)));
        print_jpeg600(fig, outFn, SAVE_DPI);
        close(fig);
    end
end

function plot_CRplusUR_composition(TCompAll, CR_BAND, UR_BANDS, phaseList, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, BandColourMap)
    if isempty(TCompAll), return; end

    comps = ["CR" UR_BANDS];

    for ph = phaseList(:).'
        phaseTag = char(ph);
        D = TCompAll(TCompAll.Phase==string(phaseTag), :);
        if isempty(D), continue; end

        ppVals = unique(D.Photoperiod_h);
        ppVals = sort(ppVals(~isnan(ppVals)));
        MeanFrac = NaN(numel(ppVals), numel(comps));

        for i = 1:numel(ppVals)
            dpp = D(D.Photoperiod_h==ppVals(i), :);
            for c = 1:numel(comps)
                cc = comps(c);
                r = dpp(dpp.Component==cc, :);
                if ~isempty(r)
                    MeanFrac(i,c) = r.MeanFrac(1);
                end
            end
        end

        fig = figure('Visible','off'); set(fig,'Position',[100 100 1400 560]);
        ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

        X = categorical(string(ppVals));
        bb = bar(X, MeanFrac, 'stacked', 'FaceColor','flat');

        for c = 1:numel(comps)
            cc = comps(c);
            if cc == "CR"
                key = char(CR_BAND);
            else
                key = char(cc);
            end
            col = [0 0 0];
            if isKey(BandColourMap, key), col = BandColourMap(key); end
            bb(c).CData = repmat(col, numel(ppVals), 1);
        end

        ylim([0 1]);
        yl = ylabel('Fraction of (CR + UR) (linear)');
        xl = xlabel('Photoperiod (h)');
        ttl = title(sprintf('Composition (CR + UR) | Phase: %s', phaseTag), 'Interpreter','none');
        set_bold_labels(ttl, xl, yl);

        lgLabels = cellstr(comps);
        lgLabels{1} = char(CR_BAND);
        lg = legend(lgLabels, 'Location','eastoutside', 'Interpreter','none');
        set(lg,'Box','off','FontName',FONT_NAME);

        outFn = fullfile(outDir, sprintf('Composition_CRplusUR_Stacked_Phase_%s.jpg', sanitise_filename(phaseTag)));
        print_jpeg600(fig, outFn, SAVE_DPI);
        close(fig);
    end
end

%% ---------------- PCA ----------------------------------------------------

function run_pca_powers_only(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, phaseList, xlsxOut, outDir, ...
    SAVE_DPI, FONT_NAME, TICK_DIR, DO_CLUSTER, PPColourMap)

    phaseTag = "All";
    if ~any(string(phaseList)==phaseTag)
        phaseTag = string(phaseList(1));
    end

    keys = {'File','SignalID','Photoperiod_h','Phase'};
    CR = BCS(BCS.Phase==phaseTag & BCS.Source==SRC_CR & BCS.BandName==CR_BAND, [keys {'MeanBandPower_log10'}]);
    CR.Properties.VariableNames{'MeanBandPower_log10'} = 'CR_Log10';
    if isempty(CR), return; end

    F = CR(:, keys);
    F.CR_Log10 = CR.CR_Log10;

    for b = 1:numel(UR_BANDS)
        ub = UR_BANDS(b);
        U = BCS(BCS.Phase==phaseTag & BCS.Source==SRC_UR & BCS.BandName==ub, [keys {'MeanBandPower_log10'}]);
        U.Properties.VariableNames{'MeanBandPower_log10'} = char("UR_" + ub);
        F = outerjoin(F, U, 'Keys', keys, 'MergeKeys', true);
    end

    feat = F.Properties.VariableNames;
    feat = feat(~ismember(feat, keys));

    X = NaN(height(F), numel(feat));
    for i = 1:numel(feat)
        X(:,i) = double(F.(feat{i}));
    end

    minFiniteFeatures = max(3, floor(0.7*size(X,2)));
    [Xz, F, feat, pcaOK] = prepare_pca_matrix_for_validated_ur(X, F, feat, minFiniteFeatures, 'PCA_PowersOnly');
    if ~pcaOK, return; end

    [coeff, score, ~, ~, explained] = pca(Xz, 'Rows','complete');
    if size(score,2) < 2 || numel(explained) < 2
        warning('PCA_PowersOnly produced fewer than two PCs after filtering/imputation. Skipping PC1-PC2 plot.');
        return;
    end

    % Write PCA tables
    try
        Tload = array2table(coeff, 'VariableNames', compose('PC%d', 1:size(coeff,2)), 'RowNames', feat);
        Tscore = array2table(score, 'VariableNames', compose('PC%d', 1:size(score,2)));
        Tscore = [F(:,keys) Tscore];
        writetable(Tscore, xlsxOut, 'Sheet', 'PCA_Powers_Scores');
        writetable(resetRowNames(Tload), xlsxOut, 'Sheet', 'PCA_Powers_Loadings');
        writetable(table(explained(:), 'VariableNames', {'ExplainedVar_percent'}), xlsxOut, 'Sheet', 'PCA_Powers_Explained');
    catch
    end

    fig = figure('Visible','off'); set(fig,'Position',[100 100 1100 560]);
    ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

    pp = double(F.Photoperiod_h);

    % v6.3.1 FIX: filter isnan on ppU, not pp
    ppU = unique(pp);
    ppU = sort(ppU(~isnan(ppU)));

    hold on;
    for i = 1:numel(ppU)
        m = (pp == ppU(i));
        key = char(string(ppU(i)));
        c = [0 0 0];
        if isKey(PPColourMap, key), c = PPColourMap(key); end
        scatter(score(m,1), score(m,2), 40, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    end
    hold off;

    xl = xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    yl = ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    ttl = title(sprintf('PCA (PowersOnly) | Phase: %s', char(phaseTag)), 'Interpreter','none');
    set_bold_labels(ttl, xl, yl);

    lg = legend(cellstr(string(ppU)), 'Location','eastoutside', 'Interpreter','none');
    set(lg,'Box','off','FontName',FONT_NAME);

    outFn = fullfile(outDir, sprintf('PCA_PowersOnly_PC1_PC2_Phase_%s.jpg', sanitise_filename(char(phaseTag))));
    print_jpeg600(fig, outFn, SAVE_DPI);
    close(fig);

    if DO_CLUSTER
        try
            [cl, kBest] = cluster_pc12(score(:,1:2));
            Tassign = table(F.File, F.SignalID, F.Photoperiod_h, F.Phase, score(:,1), score(:,2), cl, ...
                'VariableNames', {'File','SignalID','Photoperiod_h','Phase','PC1','PC2','Cluster'});
            writetable(Tassign, xlsxOut, 'Sheet', 'PCA_Powers_Clusters');

            fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 600]);
            ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);
            plot_pca_clusters(score(:,1), score(:,2), cl, PPColourMap, pp, ppU);

            xl = xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
            yl = ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
            ttl = title(sprintf('PCA (PowersOnly) clusters (k=%d) | Phase: %s', kBest, char(phaseTag)), 'Interpreter','none');
            set_bold_labels(ttl, xl, yl);

            outFn = fullfile(outDir, sprintf('PCA_PowersOnly_Clusters_PC1_PC2_Phase_%s.jpg', sanitise_filename(char(phaseTag))));
            print_jpeg600(fig, outFn, SAVE_DPI);
            close(fig);
        catch
        end
    end
end

function run_pca_deltas_only(Tpair, UR_BANDS, phaseList, xlsxOut, outDir, ...
    SAVE_DPI, FONT_NAME, TICK_DIR, DO_CLUSTER, PPColourMap)

    phaseTag = "All";
    if ~any(string(phaseList)==phaseTag)
        phaseTag = string(phaseList(1));
    end

    D = Tpair(Tpair.Phase==phaseTag, :);
    if isempty(D), return; end

    keys = {'File','SignalID','Photoperiod_h','Phase'};
    base = unique(D(:, keys));
    F = base;

    for b = 1:numel(UR_BANDS)
        ub = UR_BANDS(b);
        m = (D.UR_Band==ub);
        Db = D(m, keys);
        Db.Delta = D.Delta_log10(m);
        Db.Properties.VariableNames{'Delta'} = char("Delta_" + ub);
        F = outerjoin(F, Db, 'Keys', keys, 'MergeKeys', true);
    end

    feat = F.Properties.VariableNames;
    feat = feat(~ismember(feat, keys));

    X = NaN(height(F), numel(feat));
    for i = 1:numel(feat)
        X(:,i) = double(F.(feat{i}));
    end

    minFiniteFeatures = max(2, floor(0.7*size(X,2)));
    [Xz, F, feat, pcaOK] = prepare_pca_matrix_for_validated_ur(X, F, feat, minFiniteFeatures, 'PCA_DeltasOnly');
    if ~pcaOK, return; end

    [coeff, score, ~, ~, explained] = pca(Xz, 'Rows','complete');
    if size(score,2) < 2 || numel(explained) < 2
        warning('PCA_DeltasOnly produced fewer than two PCs after filtering/imputation. Skipping PC1-PC2 plot.');
        return;
    end

    % Write PCA tables
    try
        Tload = array2table(coeff, 'VariableNames', compose('PC%d', 1:size(coeff,2)), 'RowNames', feat);
        Tscore = array2table(score, 'VariableNames', compose('PC%d', 1:size(score,2)));
        Tscore = [F(:,keys) Tscore];
        writetable(Tscore, xlsxOut, 'Sheet', 'PCA_Deltas_Scores');
        writetable(resetRowNames(Tload), xlsxOut, 'Sheet', 'PCA_Deltas_Loadings');
        writetable(table(explained(:), 'VariableNames', {'ExplainedVar_percent'}), xlsxOut, 'Sheet', 'PCA_Deltas_Explained');
    catch
    end

    fig = figure('Visible','off'); set(fig,'Position',[100 100 1100 560]);
    ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

    pp = double(F.Photoperiod_h);

    % v6.3.1 FIX: filter isnan on ppU, not pp
    ppU = unique(pp);
    ppU = sort(ppU(~isnan(ppU)));

    hold on;
    for i = 1:numel(ppU)
        m = (pp == ppU(i));
        key = char(string(ppU(i)));
        c = [0 0 0];
        if isKey(PPColourMap, key), c = PPColourMap(key); end
        scatter(score(m,1), score(m,2), 40, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    end
    hold off;

    xl = xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    yl = ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    ttl = title(sprintf('PCA (DeltasOnly) | Phase: %s', char(phaseTag)), 'Interpreter','none');
    set_bold_labels(ttl, xl, yl);

    lg = legend(cellstr(string(ppU)), 'Location','eastoutside', 'Interpreter','none');
    set(lg,'Box','off','FontName',FONT_NAME);

    outFn = fullfile(outDir, sprintf('PCA_DeltasOnly_PC1_PC2_Phase_%s.jpg', sanitise_filename(char(phaseTag))));
    print_jpeg600(fig, outFn, SAVE_DPI);
    close(fig);

    if DO_CLUSTER
        try
            [cl, kBest] = cluster_pc12(score(:,1:2));
            Tassign = table(F.File, F.SignalID, F.Photoperiod_h, F.Phase, score(:,1), score(:,2), cl, ...
                'VariableNames', {'File','SignalID','Photoperiod_h','Phase','PC1','PC2','Cluster'});
            writetable(Tassign, xlsxOut, 'Sheet', 'PCA_Deltas_Clusters');

            fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 600]);
            ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);
            plot_pca_clusters(score(:,1), score(:,2), cl, PPColourMap, pp, ppU);

            xl = xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
            yl = ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
            ttl = title(sprintf('PCA (DeltasOnly) clusters (k=%d) | Phase: %s', kBest, char(phaseTag)), 'Interpreter','none');
            set_bold_labels(ttl, xl, yl);

            outFn = fullfile(outDir, sprintf('PCA_DeltasOnly_Clusters_PC1_PC2_Phase_%s.jpg', sanitise_filename(char(phaseTag))));
            print_jpeg600(fig, outFn, SAVE_DPI);
            close(fig);
        catch
        end
    end
end


function [Xz, Fout, featOut, okFlag] = prepare_pca_matrix_for_validated_ur(X, F, feat, minFiniteFeatures, label)
    % Robust PCA pre-processing for the validated-Raw UR pipeline.
    %
    % The validated UR matrix is often sparse because not every mouse has
    % every validated UR band in every photoperiod. A strict complete-case
    % PCA can therefore collapse to one PC or too few rows. This helper:
    %   1) removes unusable features;
    %   2) keeps rows with a minimum number of observed features;
    %   3) z-scores using available values;
    %   4) mean-imputes remaining missing values after z-scoring, i.e. NaN=0.
    %
    % This preserves PCA as an exploratory visualisation while preventing
    % the PC1-PC2 plotting step from crashing when complete cases are sparse.

    okFlag = false;
    Xz = [];
    Fout = F;
    featOut = feat;

    if isempty(X) || isempty(feat)
        warning('%s skipped: empty PCA matrix.', label);
        return;
    end

    finitePerCol = sum(isfinite(X), 1);
    colStd = nan(1, size(X,2));
    for j = 1:size(X,2)
        vals = X(isfinite(X(:,j)), j);
        if numel(vals) >= 2
            colStd(j) = std(vals, 0);
        end
    end

    keepCol = finitePerCol >= 3 & isfinite(colStd) & colStd > 0;
    X = X(:, keepCol);
    featOut = feat(keepCol);

    if size(X,2) < 2
        warning('%s skipped: fewer than two usable PCA features after filtering.', label);
        return;
    end

    minFiniteFeatures = min(minFiniteFeatures, size(X,2));
    keepRow = sum(isfinite(X), 2) >= minFiniteFeatures;
    X = X(keepRow, :);
    Fout = F(keepRow, :);

    if size(X,1) < 3
        warning('%s skipped: fewer than three usable rows after PCA filtering.', label);
        return;
    end

    mu = nan(1, size(X,2));
    sig = nan(1, size(X,2));

    for j = 1:size(X,2)
        vals = X(isfinite(X(:,j)), j);
        mu(j) = mean(vals, 'omitnan');
        sig(j) = std(vals, 0, 'omitnan');
        if ~isfinite(sig(j)) || sig(j) == 0
            sig(j) = 1;
        end
    end

    Xz = (X - mu) ./ sig;
    Xz(~isfinite(Xz)) = 0;  % mean imputation after z-scoring

    if size(Xz,1) < 3 || size(Xz,2) < 2
        warning('%s skipped: PCA matrix too small after standardisation.', label);
        return;
    end

    okFlag = true;
end

function [cl, kBest] = cluster_pc12(X2)
    % Simple clustering on PC1-PC2 for visualisation.
    % Try k=2..6 and choose best silhouette. If silhouette fails, fallback k=2.
    kList = 2:min(6, max(2, size(X2,1)-1));
    bestS = -Inf; kBest = 2; cl = ones(size(X2,1),1);

    for k = kList
        try
            c = kmeans(X2, k, 'Replicates', 10, 'MaxIter', 500, 'Display', 'off');
            s = mean(silhouette(X2, c), 'omitnan');
            if s > bestS
                bestS = s;
                kBest = k;
                cl = c;
            end
        catch
        end
    end
end

function plot_pca_clusters(pc1, pc2, cl, PPColourMap, pp, ppU)
    % Plot points coloured by photoperiod, with cluster hulls.
    hold on;

    for i = 1:numel(ppU)
        m = (pp == ppU(i));
        key = char(string(ppU(i)));
        c = [0 0 0];
        if isKey(PPColourMap, key), c = PPColourMap(key); end
        scatter(pc1(m), pc2(m), 40, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    end

    u = unique(cl);
    for j = 1:numel(u)
        m = (cl == u(j));
        if sum(m) >= 4
            try
                K = convhull(pc1(m), pc2(m));
                plot(pc1(m(K)), pc2(m(K)), 'k-', 'LineWidth', 1.2);
            catch
            end
        end
    end

    hold off;
end

function T = resetRowNames(T)
    if isempty(T.Properties.RowNames), return; end
    rn = T.Properties.RowNames;
    T = addvars(T, string(rn), 'Before', 1, 'NewVariableNames', 'Feature');
    T.Properties.RowNames = {};
end

%% ---------------- TS: TRUE COVERAGE + DOMINANCE --------------------------

function Tcov = compute_true_coverage(BP, RP, RPOW, phasesToUse)
    phasesToUse = string(phasesToUse);
    rows = {};
    hdr = {'SeriesType','File','SignalID','Photoperiod_h','Phase','BandName','Source','Coverage','N_valid','N_total'};

    rows = [rows; truecov_rows_one(BP, phasesToUse, 'BandPower')]; %#ok<AGROW>
    rows = [rows; truecov_rows_one(RP, phasesToUse, 'RidgePeriod')]; %#ok<AGROW>
    rows = [rows; truecov_rows_one(RPOW, phasesToUse, 'RidgePower')]; %#ok<AGROW>

    if isempty(rows)
        Tcov = table();
    else
        Tcov = cell2table(rows, 'VariableNames', hdr);
        Tcov.SeriesType = string(Tcov.SeriesType);
        Tcov.File = string(Tcov.File);
        Tcov.SignalID = string(Tcov.SignalID);
        Tcov.Phase = string(Tcov.Phase);
        Tcov.BandName = string(Tcov.BandName);
        Tcov.Source = string(Tcov.Source);
    end
end

function rows = truecov_rows_one(T, phasesToUse, seriesType)
    rows = {};
    if isempty(T), return; end
    req = ["File","SignalID","Source","BandName","Phase","Value","ValidFlag"];
    if any(~ismember(req, string(T.Properties.VariableNames))), return; end

    T = T(ismember(T.Phase, phasesToUse), :);
    if isempty(T), return; end

    % total points per series
    pp = parse_photoperiod_from_file(T.File);
    G = findgroups(T.File, T.SignalID, pp, T.Phase, T.BandName, T.Source);

    nTot = splitapply(@(v) numel(v), T.Value, G);

    isValid = (double(T.ValidFlag)==1) & isfinite(T.Value);
    nVal = splitapply(@(m) sum(m), isValid, G);

    cov = nVal ./ max(nTot,1);

    file = splitapply(@(x) x(1), T.File, G);
    sid  = splitapply(@(x) x(1), T.SignalID, G);
    pp1  = splitapply(@(x) x(1), pp, G);
    ph   = splitapply(@(x) x(1), T.Phase, G);
    bn   = splitapply(@(x) x(1), T.BandName, G);
    src  = splitapply(@(x) x(1), T.Source, G);

    for i = 1:numel(cov)
        rows(end+1,:) = {seriesType, file(i), sid(i), pp1(i), ph(i), bn(i), src(i), cov(i), nVal(i), nTot(i)}; %#ok<AGROW>
    end
end

function Tdom = compute_dominance_occupancy_BP_truecov(BP, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, phasesToUse)
    phasesToUse = string(phasesToUse);
    req = ["File","SignalID","Source","BandName","Phase","Time_days","Value","ValidFlag"];
    if any(~ismember(req, string(BP.Properties.VariableNames))), Tdom = table(); return; end

    BP = BP(ismember(BP.Phase, phasesToUse), :);

    CR = BP(BP.Source==SRC_CR & BP.BandName==CR_BAND, :);
    CR = CR(:, {'File','SignalID','Phase','Time_days','Value','ValidFlag'});
    CR.Properties.VariableNames{'Value'} = 'CR_Value';
    CR.Properties.VariableNames{'ValidFlag'} = 'CR_OK';

    rows = {};
    hdr = {'File','SignalID','Photoperiod_h','Phase','UR_Band','Occ_URgtCR','N_validPair','N_totalPair','Coverage_CR','Coverage_UR'};

    for b = 1:numel(UR_BANDS)
        ub = UR_BANDS(b);
        U = BP(BP.Source==SRC_UR & BP.BandName==ub, :);
        if isempty(U), continue; end
        U = U(:, {'File','SignalID','Phase','Time_days','Value','ValidFlag'});
        U.Properties.VariableNames{'Value'} = 'UR_Value';
        U.Properties.VariableNames{'ValidFlag'} = 'UR_OK';

        J = innerjoin(U, CR, 'Keys', {'File','SignalID','Phase','Time_days'});
        if isempty(J), continue; end

        nTot = height(J);

        okPair = (double(J.UR_OK)==1) & (double(J.CR_OK)==1) & isfinite(J.UR_Value) & isfinite(J.CR_Value);
        nVal = sum(okPair);
        if nVal == 0, continue; end

        occ = mean(double(J.UR_Value(okPair) > J.CR_Value(okPair)), 'omitnan');

        % Coverage within paired grid
        covCR = mean((double(J.CR_OK)==1) & isfinite(J.CR_Value));
        covUR = mean((double(J.UR_OK)==1) & isfinite(J.UR_Value));

        pp = parse_photoperiod_from_file(J.File);
        if isempty(pp), pp = NaN; else, pp = pp(1); end

        rows(end+1,:) = {J.File(1), J.SignalID(1), pp, J.Phase(1), ub, occ, nVal, nTot, covCR, covUR}; %#ok<AGROW>
    end

    if isempty(rows)
        Tdom = table();
    else
        Tdom = cell2table(rows, 'VariableNames', hdr);
        Tdom.File = string(Tdom.File);
        Tdom.SignalID = string(Tdom.SignalID);
        Tdom.Phase = string(Tdom.Phase);
        Tdom.UR_Band = string(Tdom.UR_Band);
    end
end

function Tstab = compute_stability_from_long_truecov(Tlong, phasesToUse, seriesType)
    phasesToUse = string(phasesToUse);
    req = ["File","SignalID","Source","BandName","Phase","Value","ValidFlag"];
    if any(~ismember(req, string(Tlong.Properties.VariableNames))), Tstab = table(); return; end

    Tlong = Tlong(ismember(Tlong.Phase, phasesToUse), :);
    if isempty(Tlong), Tstab = table(); return; end

    pp = parse_photoperiod_from_file(Tlong.File);

    % totals per series
    G = findgroups(Tlong.File, Tlong.SignalID, pp, Tlong.Phase, Tlong.BandName, Tlong.Source);
    nTot = splitapply(@(v) numel(v), Tlong.Value, G);

    isValid = (double(Tlong.ValidFlag)==1) & isfinite(Tlong.Value);

    sdV  = splitapply(@(v,m) std(v(m),'omitnan'),  Tlong.Value, isValid, G);
    madV = splitapply(@(v,m) mad(v(m),1),         Tlong.Value, isValid, G);
    nVal = splitapply(@(m) sum(m), isValid, G);

    cov = nVal ./ max(nTot,1);

    file = splitapply(@(x) x(1), Tlong.File, G);
    sid  = splitapply(@(x) x(1), Tlong.SignalID, G);
    pp1  = splitapply(@(x) x(1), pp, G);
    ph   = splitapply(@(x) x(1), Tlong.Phase, G);
    bn   = splitapply(@(x) x(1), Tlong.BandName, G);
    src  = splitapply(@(x) x(1), Tlong.Source, G);

    Tstab = table( ...
        repmat(string(seriesType), numel(sdV), 1), ...
        string(file), string(sid), pp1, string(ph), string(bn), string(src), ...
        sdV, madV, cov, nVal, nTot, ...
        'VariableNames', {'SeriesType','File','SignalID','Photoperiod_h','Phase','BandName','Source', ...
                          'SD_overTime','MAD_overTime','Coverage','N_valid','N_total'});
end

function pp = parse_photoperiod_from_file(fileStr)
    s = string(fileStr);
    pp = NaN(numel(s),1);
    for i = 1:numel(s)
        tok = regexp(char(s(i)), 'L(\d{1,2})', 'tokens', 'once');
        if ~isempty(tok)
            pp(i) = str2double(tok{1});
        end
    end
end

function plot_coverage_fraction_bars(Tcov, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, PPColourMap)
    if isempty(Tcov), return; end
    phases = unique(Tcov.Phase);

    for ph = phases(:).'
        D = Tcov(Tcov.Phase==ph & Tcov.SeriesType=="BandPower", :);
        if isempty(D), continue; end

        bands = unique(D.BandName);
        for bn = bands(:).'
            Db = D(D.BandName==bn, :);
            if isempty(Db), continue; end

            G = findgroups(Db.Photoperiod_h);
            pp = splitapply(@(x) x(1), Db.Photoperiod_h, G);
            m  = splitapply(@(x) mean(x,'omitnan'), Db.Coverage, G);
            s  = splitapply(@(x) std(x,'omitnan'),  Db.Coverage, G);

            [pp, ix] = sort(pp);
            m = m(ix); s = s(ix);

            fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
            ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

            X = categorical(string(pp));
            bb = bar(X, m, 'FaceColor','flat');
            for i = 1:numel(pp)
                key = char(string(pp(i)));
                if isKey(PPColourMap, key)
                    bb.CData(i,:) = PPColourMap(key);
                end
            end

            hold on;
            errorbar(bb.XEndPoints, m, s, 'k', 'LineStyle','none', 'LineWidth', 1.2);
            hold off;

            ylim([0 1]);
            yl = ylabel('Coverage fraction (ValidFlag & finite)');
            xl = xlabel('Photoperiod (h)');
            ttl = title(sprintf('Coverage: %s | Phase: %s', char(bn), char(ph)), 'Interpreter','none');
            set_bold_labels(ttl, xl, yl);

            outFn = fullfile(outDir, sprintf('Coverage_%s_Phase_%s.jpg', sanitise_filename(char(bn)), sanitise_filename(char(ph))));
            print_jpeg600(fig, outFn, SAVE_DPI);
            close(fig);
        end
    end
end

function plot_dominance_bars(Tdom, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, UR_BANDS, PPColourMap)
    if isempty(Tdom), return; end
    phases = unique(Tdom.Phase);

    for ph = phases(:).'
        D = Tdom(Tdom.Phase==ph, :);
        if isempty(D), continue; end

        for b = 1:numel(UR_BANDS)
            ub = UR_BANDS(b);
            Db = D(D.UR_Band==ub, :);
            if isempty(Db), continue; end

            G = findgroups(Db.Photoperiod_h);
            pp = splitapply(@(x) x(1), Db.Photoperiod_h, G);
            m  = splitapply(@(x) mean(x,'omitnan'), Db.Occ_URgtCR, G);
            s  = splitapply(@(x) std(x,'omitnan'),  Db.Occ_URgtCR, G);

            [pp, ix] = sort(pp);
            m = m(ix); s = s(ix);

            fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
            ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

            X = categorical(string(pp));
            bb = bar(X, m, 'FaceColor','flat');
            for i = 1:numel(pp)
                key = char(string(pp(i)));
                if isKey(PPColourMap, key)
                    bb.CData(i,:) = PPColourMap(key);
                end
            end

            hold on;
            errorbar(bb.XEndPoints, m, s, 'k', 'LineStyle','none', 'LineWidth', 1.2);
            yline(0.5,'k:','LineWidth',1.0);
            hold off;

            yl = ylabel('P(UR_{band} > CR) (BandPower)');
            xl = xlabel('Photoperiod (h)');
            ttl = title(sprintf('Dominance: %s | Phase: %s', char(ub), char(ph)), 'Interpreter','none');
            set_bold_labels(ttl, xl, yl);

            outFn = fullfile(outDir, sprintf('Dominance_%s_Phase_%s.jpg', sanitise_filename(char(ub)), sanitise_filename(char(ph))));
            print_jpeg600(fig, outFn, SAVE_DPI);
            close(fig);
        end
    end
end

function plot_stability_bars(Tstab, outDir, SAVE_DPI, FONT_NAME, TICK_DIR, labelStr, PPColourMap)
    if isempty(Tstab), return; end
    phases = unique(Tstab.Phase);

    for ph = phases(:).'
        D = Tstab(Tstab.Phase==ph, :);
        if isempty(D), continue; end

        bands = unique(D.BandName);
        for bn = bands(:).'
            Db = D(D.BandName==bn, :);
            if isempty(Db), continue; end

            G = findgroups(Db.Photoperiod_h);
            pp = splitapply(@(x) x(1), Db.Photoperiod_h, G);

            mSD  = splitapply(@(x) mean(x,'omitnan'), Db.SD_overTime, G);
            sSD  = splitapply(@(x) std(x,'omitnan'),  Db.SD_overTime, G);
            mMAD = splitapply(@(x) mean(x,'omitnan'), Db.MAD_overTime, G);
            sMAD = splitapply(@(x) std(x,'omitnan'),  Db.MAD_overTime, G);

            [pp, ix] = sort(pp);
            mSD = mSD(ix); sSD = sSD(ix);
            mMAD = mMAD(ix); sMAD = sMAD(ix);

            fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
            ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

            X = categorical(string(pp));
            bb = bar(X, mSD, 'FaceColor','flat');
            for i = 1:numel(pp)
                key = char(string(pp(i)));
                if isKey(PPColourMap, key)
                    bb.CData(i,:) = PPColourMap(key);
                end
            end

            hold on; errorbar(bb.XEndPoints, mSD, sSD, 'k', 'LineStyle','none', 'LineWidth', 1.2); hold off;

            yl = ylabel('SD over time');
            xl = xlabel('Photoperiod (h)');
            ttl = title(sprintf('%s stability (SD): %s | Phase: %s', labelStr, char(bn), char(ph)), 'Interpreter','none');
            set_bold_labels(ttl, xl, yl);

            outFn = fullfile(outDir, sprintf('%s_SD_%s_Phase_%s.jpg', sanitise_filename(labelStr), sanitise_filename(char(bn)), sanitise_filename(char(ph))));
            print_jpeg600(fig, outFn, SAVE_DPI);
            close(fig);

            fig = figure('Visible','off'); set(fig,'Position',[100 100 1200 520]);
            ax = gca; apply_axes_style(ax, FONT_NAME, TICK_DIR);

            X = categorical(string(pp));
            bb = bar(X, mMAD, 'FaceColor','flat');
            for i = 1:numel(pp)
                key = char(string(pp(i)));
                if isKey(PPColourMap, key)
                    bb.CData(i,:) = PPColourMap(key);
                end
            end

            hold on; errorbar(bb.XEndPoints, mMAD, sMAD, 'k', 'LineStyle','none', 'LineWidth', 1.2); hold off;

            yl = ylabel('MAD over time');
            xl = xlabel('Photoperiod (h)');
            ttl = title(sprintf('%s stability (MAD): %s | Phase: %s', labelStr, char(bn), char(ph)), 'Interpreter','none');
            set_bold_labels(ttl, xl, yl);

            outFn = fullfile(outDir, sprintf('%s_MAD_%s_Phase_%s.jpg', sanitise_filename(labelStr), sanitise_filename(char(bn)), sanitise_filename(char(ph))));
            print_jpeg600(fig, outFn, SAVE_DPI);
            close(fig);
        end
    end
end

%% ========================================================================
% v7 validated-Raw helper functions
% ========================================================================

function C = harmonise_carryforward_table(C)
    if isempty(C), return; end
    vars = C.Properties.VariableNames;
    strVars = {'File','SignalID','ConditionParsed','Phase','BandName','RawCandidateID','PrimaryHSubMode', ...
               'PrimaryHSubCandidateID','MatchStatus','QCFlag','FullLadderSensitivityStatus','FinalValidationClass'};
    for i = 1:numel(strVars)
        if ismember(strVars{i}, vars)
            C.(strVars{i}) = string(C.(strVars{i}));
        end
    end
    if ismember('Photoperiod_h', vars)
        C.Photoperiod_h = double(C.Photoperiod_h);
    end
    if ismember('CarryForward', vars)
        C.CarryForward = logicalish(C.CarryForward);
        C = C(C.CarryForward, :);
    end
end

function T = annotate_BCS_with_validation(T, CarryForward, UR_BANDS, SRC_UR, applyAllPhase)
    if isempty(T), return; end
    T.Source   = string(T.Source);
    T.BandName = string(T.BandName);
    T.File     = string(T.File);
    T.SignalID = string(T.SignalID);
    T.Phase    = string(T.Phase);
    if ismember('Photoperiod_h', T.Properties.VariableNames)
        T.Photoperiod_h = double(T.Photoperiod_h);
    end

    n = height(T);
    T.ValidatedRawUR = false(n,1);
    T.ValidationStatus = repmat("NotApplicable", n, 1);
    T.ValidationSourcePhase = strings(n,1); T.ValidationSourcePhase(:) = missing;
    T.ValidationRawCandidateID = strings(n,1); T.ValidationRawCandidateID(:) = missing;
    T.ValidationFullLadderSensitivity = strings(n,1); T.ValidationFullLadderSensitivity(:) = missing;
    T.ValidationFinalClass = strings(n,1); T.ValidationFinalClass(:) = missing;

    isRawUR = T.Source==SRC_UR & ismember(T.BandName, UR_BANDS);
    T.ValidationStatus(isRawUR) = "RawUR_NotValidated";

    if isempty(CarryForward)
        return;
    end

    % Primary validation is normally done on Phase="All". By default, this
    % all-phase decision is propagated to Light/Dark summaries for the same
    % File x SignalID x Photoperiod x Band candidate.
    keysCarryAll = make_validation_keys(CarryForward, false);
    keysCarryExact = make_validation_keys(CarryForward, true);

    for r = find(isRawUR).'
        keyAll = make_one_validation_key(T.File(r), T.SignalID(r), T.Photoperiod_h(r), T.BandName(r), "");
        idx = find(keysCarryAll == keyAll, 1, 'first');

        sourcePhase = "All";
        if ~applyAllPhase || isempty(idx)
            keyExact = make_one_validation_key(T.File(r), T.SignalID(r), T.Photoperiod_h(r), T.BandName(r), T.Phase(r));
            idx = find(keysCarryExact == keyExact, 1, 'first');
            sourcePhase = T.Phase(r);
        end

        if ~isempty(idx)
            T.ValidatedRawUR(r) = true;
            T.ValidationStatus(r) = "CarryForward_ValidatedRawUR";
            T.ValidationSourcePhase(r) = sourcePhase;
            if ismember('RawCandidateID', CarryForward.Properties.VariableNames)
                T.ValidationRawCandidateID(r) = string(CarryForward.RawCandidateID(idx));
            end
            if ismember('FullLadderSensitivityStatus', CarryForward.Properties.VariableNames)
                T.ValidationFullLadderSensitivity(r) = string(CarryForward.FullLadderSensitivityStatus(idx));
            end
            if ismember('FinalValidationClass', CarryForward.Properties.VariableNames)
                T.ValidationFinalClass(r) = string(CarryForward.FinalValidationClass(idx));
            end
        end
    end
end

function Tused = filter_BCS_for_validated_raw_analysis(T, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, requireValidatedUR)
    if isempty(T), Tused = T; return; end
    T.Source = string(T.Source); T.BandName = string(T.BandName);
    keepCR = T.Source==SRC_CR & T.BandName==CR_BAND;
    keepUR = T.Source==SRC_UR & ismember(T.BandName, UR_BANDS);
    if requireValidatedUR && ismember('ValidatedRawUR', T.Properties.VariableNames)
        keepUR = keepUR & T.ValidatedRawUR;
    end
    Tused = T(keepCR | keepUR, :);
end


function T = filter_time_days_minimum(T, minDay)
    if isempty(T) || ~istable(T) || ~ismember('Time_days', T.Properties.VariableNames)
        return;
    end
    try
        T.Time_days = double(T.Time_days);
        T = T(isfinite(T.Time_days) & T.Time_days >= minDay, :);
    catch
        % Leave table unchanged if conversion fails.
    end
end

function T = filter_long_for_validated_raw_analysis(T, CarryForward, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, useValidation, applyAllPhase)
    if isempty(T), return; end
    req = {'File','SignalID','Source','BandName','Phase'};
    if ~all(ismember(req, T.Properties.VariableNames))
        return;
    end

    T.File     = string(T.File);
    T.SignalID = string(T.SignalID);
    T.Source   = string(T.Source);
    T.BandName = string(T.BandName);
    T.Phase    = string(T.Phase);

    keepCR = T.Source==SRC_CR & T.BandName==CR_BAND;
    keepUR = T.Source==SRC_UR & ismember(T.BandName, UR_BANDS);

    if useValidation && ~isempty(CarryForward)
        % Legacy TS tables lack Photoperiod_h and CandidateID, so the validated
        % gate is applied by File x SignalID x Band. If validation was done on
        % Phase=All, that validated decision is propagated to Light/Dark TS rows.
        keysCarryNoPP = make_validation_keys_no_pp(CarryForward, false);
        keysTNoPP = make_validation_keys_no_pp(T, false);
        validUR = false(height(T),1);
        for r = find(keepUR).'
            if any(keysCarryNoPP == keysTNoPP(r))
                validUR(r) = true;
            elseif ~applyAllPhase
                keysCarryExact = make_validation_keys_no_pp(CarryForward, true);
                keysTExact = make_validation_keys_no_pp(T, true);
                validUR(r) = any(keysCarryExact == keysTExact(r));
            end
        end
        keepUR = keepUR & validUR;
    end

    T = T(keepCR | keepUR, :);
end

function keys = make_validation_keys(T, includePhase)
    if isempty(T), keys = strings(0,1); return; end
    n = height(T); keys = strings(n,1);
    for i = 1:n
        phase = "";
        if includePhase && ismember('Phase', T.Properties.VariableNames)
            phase = string(T.Phase(i));
        end
        pp = NaN;
        if ismember('Photoperiod_h', T.Properties.VariableNames), pp = double(T.Photoperiod_h(i)); end
        keys(i) = make_one_validation_key(string(T.File(i)), string(T.SignalID(i)), pp, string(T.BandName(i)), phase);
    end
end

function keys = make_validation_keys_no_pp(T, includePhase)
    if isempty(T), keys = strings(0,1); return; end
    n = height(T); keys = strings(n,1);
    for i = 1:n
        phase = "";
        if includePhase && ismember('Phase', T.Properties.VariableNames)
            phase = string(T.Phase(i));
        end
        keys(i) = lower(strjoin([string(T.File(i)), string(T.SignalID(i)), string(T.BandName(i)), phase], "|"));
    end
end

function key = make_one_validation_key(file, signalID, photoperiod_h, bandName, phase)
    if isempty(phase), phase = ""; end
    key = lower(strjoin([string(file), string(signalID), sprintf('%.10g', double(photoperiod_h)), string(bandName), string(phase)], "|"));
end

function x = logicalish(x)
    if islogical(x)
        return;
    elseif isnumeric(x)
        x = x(:) ~= 0;
    else
        s = lower(strtrim(string(x(:))));
        x = ismember(s, ["true","1","yes","y","pass","passed"]);
    end
end
