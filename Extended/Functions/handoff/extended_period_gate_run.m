function validationMap = extended_period_gate_run(handoffDir, cfg)
%EXTENDED_PERIOD_GATE_RUN CarryForward gate: Raw UR vs SEL_P360 (±15%).
%
%   validationMap = extended_period_gate_run()
%   validationMap = extended_period_gate_run(handoffDir)
%   validationMap = extended_period_gate_run(handoffDir, cfg)
%
%   Kent C logic modularised for Extended E2. Raw remains the biological
%   signal; Selective-HSub (PRIMARY_HSUB_MODE, default SEL_P360) is validation
%   only. Parallel to — not interchangeable with — Core Script 3 residual-CWT UR.
%
%   Inputs
%     handoffDir  AcrossPhotoperiod_Input folder with WP_Summary__*.mat.
%                 If omitted / empty, uigetdir (legacy interactive behaviour).
%     cfg         Optional struct from extended_defaults(); missing fields filled.
%
%   Output
%     validationMap  Struct written to HSubSupported_PeriodMap.mat (also returned).
%                    Empty if user cancels folder selection.
%
%   See also: run_extended_period_gate, extended_defaults, HANDOFF_SCHEMA_EXTENDED.md

    SCRIPT_NAME    = 'extended_period_gate_run';
    SCRIPT_VERSION = '2.0-E2';
    RUN_TIMESTAMP  = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    if nargin < 2 || isempty(cfg)
        if exist('extended_defaults', 'file') == 2
            cfg = extended_defaults();
        else
            cfg = struct();
        end
    end
    cfg = fill_gate_cfg_(cfg);

    if nargin < 1 || isempty(handoffDir)
        fprintf('\n%s\n', SCRIPT_NAME);
        fprintf('Select the AcrossPhotoperiod_Input folder containing WP_Summary__*.mat files.\n');
        handoffDir = uigetdir(pwd, ...
            'Select AcrossPhotoperiod_Input folder containing WP_Summary__*.mat');
        if isequal(handoffDir, 0)
            fprintf('No folder selected. Exiting.\n');
            validationMap = [];
            return;
        end
    end
    handoffDir = char(handoffDir);

    PRIMARY_HSUB_MODE      = string(cfg.hsub.primaryMode);
    SECONDARY_HSUB_MODES   = string(cfg.hsub.secondaryModes(:));
    SENSITIVITY_HSUB_MODES = string(cfg.hsub.sensitivityModes(:));
    PRIMARY_PHASE          = string(cfg.carryForward.primaryPhase);
    UR_BANDS               = resolve_ur_band_names_(cfg);
    HARMONIC_SENSITIVE_BANDS = string(cfg.bands.harmonicSensitive(:));

    PERIOD_TOLERANCE_FRAC = cfg.carryForward.periodToleranceFrac;
    PERIOD_TOLERANCE_LOG2 = log2(1 + PERIOD_TOLERANCE_FRAC);
    MIN_RAW_RIDGE_COVERAGE  = cfg.carryForward.minRawRidgeCoverage;
    MIN_RAW_COI_VALID_FRAC  = cfg.carryForward.minRawCOIValidFrac;
    MIN_HSUB_RIDGE_COVERAGE = cfg.carryForward.minHSubRidgeCoverage;
    MIN_HSUB_COI_VALID_FRAC = cfg.carryForward.minHSubCOIValidFrac;
    REQUIRE_RAW_PASSQC  = cfg.carryForward.requireRawPassQC;
    REQUIRE_HSUB_PASSQC = cfg.carryForward.requireHSubPassQC;
    SAVE_DPI = cfg.plot.saveDpi;
    FIG_EXT  = cfg.plot.figExt;

    sumFiles = dir(fullfile(handoffDir, 'WP_Summary__*.mat'));
    if isempty(sumFiles)
        error('extended_period_gate_run:NoSummaries', ...
            'No WP_Summary__*.mat files found in: %s', handoffDir);
    end

    outRoot = fullfile(handoffDir, 'RawVsSelectiveHSub_PeriodValidation');
    figDir  = fullfile(outRoot, 'Figures');
    logDir  = fullfile(outRoot, 'Logs');
    extended_period_gate_ensure_dir(outRoot);
    extended_period_gate_ensure_dir(figDir);
    extended_period_gate_ensure_dir(logDir);

    logPath = fullfile(logDir, sprintf('%s_Log_%s.txt', SCRIPT_NAME, datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupObj = onCleanup(@() extended_period_gate_fclose_if_open(LOG)); %#ok<NASGU>
    extended_period_gate_log(LOG, '%s started at %s', SCRIPT_NAME, RUN_TIMESTAMP);
    extended_period_gate_log(LOG, 'Input folder: %s', handoffDir);
    extended_period_gate_log(LOG, 'Found %d WP_Summary files.', numel(sumFiles));
    extended_period_gate_log(LOG, 'PrimaryHSubMode=%s | PeriodToleranceFrac=%.4f', ...
        PRIMARY_HSUB_MODE, PERIOD_TOLERANCE_FRAC);

    [allCandidates, LoadSummary] = extended_period_gate_load(sumFiles, LOG);
    if isempty(allCandidates) || height(allCandidates) == 0
        error('extended_period_gate_run:NoCandidates', ...
            'No PeriodCandidates_Long rows were loaded. See log: %s', logPath);
    end

    allCandidates.GlobalCandidateRow = (1:height(allCandidates))';
    allCandidates = movevars(allCandidates, 'GlobalCandidateRow', 'Before', 1);

    T = allCandidates;
    T.SourceNorm = upper(string(T.Source));
    T.ModeNorm   = upper(string(T.HSubResidualMode));
    T.PhaseNorm  = string(T.Phase);
    T.BandNorm   = string(T.BandName);

    isUR = ismember(T.BandNorm, UR_BANDS);
    isPhase = strcmpi(T.PhaseNorm, PRIMARY_PHASE);
    Tval = T(isUR & isPhase, :);

    if isempty(Tval) || height(Tval) == 0
        error('extended_period_gate_run:NoURPhase', ...
            'No ultradian candidate rows found for Phase=%s.', PRIMARY_PHASE);
    end

    rawMask = strcmpi(Tval.SourceNorm, 'RAW');
    hsubMask = strcmpi(Tval.SourceNorm, 'HSUB');
    primaryMask = hsubMask & strcmpi(Tval.ModeNorm, upper(PRIMARY_HSUB_MODE));

    RawCandidates = Tval(rawMask, :);
    PrimaryHSubCandidates = Tval(primaryMask, :);

    if isempty(RawCandidates) || height(RawCandidates) == 0
        error('extended_period_gate_run:NoRaw', ...
            'No Raw ultradian candidates found in PeriodCandidates_Long.');
    end
    if isempty(PrimaryHSubCandidates) || height(PrimaryHSubCandidates) == 0
        warning('extended_period_gate_run:NoPrimaryHSub', ...
            'No primary HSub candidates found for mode %s. No candidates will be carried forward.', ...
            PRIMARY_HSUB_MODE);
    end

    RawCandidates.EligibleRaw = extended_period_gate_passes_qc( ...
        RawCandidates, MIN_RAW_RIDGE_COVERAGE, MIN_RAW_COI_VALID_FRAC, REQUIRE_RAW_PASSQC);
    PrimaryHSubCandidates.EligibleHSub = extended_period_gate_passes_qc( ...
        PrimaryHSubCandidates, MIN_HSUB_RIDGE_COVERAGE, MIN_HSUB_COI_VALID_FRAC, REQUIRE_HSUB_PASSQC);

    extended_period_gate_log(LOG, 'Raw UR candidate rows: %d', height(RawCandidates));
    extended_period_gate_log(LOG, 'Primary HSub UR candidate rows (%s): %d', ...
        PRIMARY_HSUB_MODE, height(PrimaryHSubCandidates));
    extended_period_gate_log(LOG, 'Eligible Raw candidates: %d', sum(RawCandidates.EligibleRaw));
    extended_period_gate_log(LOG, 'Eligible primary HSub candidates: %d', ...
        sum(PrimaryHSubCandidates.EligibleHSub));

    Matched = extended_period_gate_match(RawCandidates, PrimaryHSubCandidates, ...
        PRIMARY_HSUB_MODE, PERIOD_TOLERANCE_LOG2, PERIOD_TOLERANCE_FRAC, ...
        HARMONIC_SENSITIVE_BANDS);

    modeList = unique([SECONDARY_HSUB_MODES(:); SENSITIVITY_HSUB_MODES(:)], 'stable');
    for m = 1:numel(modeList)
        mode = modeList(m);
        Hm = Tval(hsubMask & strcmpi(Tval.ModeNorm, upper(mode)), :);
        if ~isempty(Hm)
            Hm.EligibleHSub = extended_period_gate_passes_qc(Hm, ...
                MIN_HSUB_RIDGE_COVERAGE, MIN_HSUB_COI_VALID_FRAC, REQUIRE_HSUB_PASSQC);
        end
        Matched = extended_period_gate_add_mode_columns(Matched, Hm, mode, PERIOD_TOLERANCE_LOG2);
    end

    Matched = extended_period_gate_add_fl_sensitivity(Matched, SENSITIVITY_HSUB_MODES);
    Matched = extended_period_gate_add_final_class(Matched);

    CarryForward = Matched(Matched.CarryForward, :);
    RawOnly_NotCarriedForward = Matched(~Matched.CarryForward, :);

    HSubOnly_ResidualFeatures = extended_period_gate_hsub_only( ...
        RawCandidates, PrimaryHSubCandidates, PERIOD_TOLERANCE_LOG2, PRIMARY_HSUB_MODE);
    FullLadder_Sensitivity = extended_period_gate_fl_sensitivity_table(Matched, SENSITIVITY_HSUB_MODES);
    Retention_ByBand = extended_period_gate_retention(Matched, "BandName");
    Retention_ByPhotoperiod = extended_period_gate_retention(Matched, "Photoperiod_h");
    Retention_ByPhotoperiodBand = extended_period_gate_retention(Matched, ["Photoperiod_h", "BandName"]);
    QC_Flags = extended_period_gate_qc_flags(Matched);

    extended_period_gate_log(LOG, 'Carry-forward candidates: %d', height(CarryForward));
    extended_period_gate_log(LOG, 'Raw-only/not-carried-forward candidates: %d', ...
        height(RawOnly_NotCarriedForward));
    extended_period_gate_log(LOG, 'HSub-only residual features: %d', ...
        height(HSubOnly_ResidualFeatures));

    Settings = table();
    Settings.ScriptName = string(SCRIPT_NAME);
    Settings.ScriptVersion = string(SCRIPT_VERSION);
    Settings.Timestamp = RUN_TIMESTAMP;
    Settings.InputFolder = string(handoffDir);
    Settings.PrimaryHSubMode = PRIMARY_HSUB_MODE;
    Settings.SecondaryHSubModes = join(SECONDARY_HSUB_MODES, '; ');
    Settings.SensitivityHSubModes = join(SENSITIVITY_HSUB_MODES, '; ');
    Settings.PrimaryPhase = PRIMARY_PHASE;
    Settings.PeriodToleranceFrac = PERIOD_TOLERANCE_FRAC;
    Settings.PeriodToleranceLog2 = PERIOD_TOLERANCE_LOG2;
    Settings.MinRawRidgeCoverage = MIN_RAW_RIDGE_COVERAGE;
    Settings.MinRawCOIValidFrac = MIN_RAW_COI_VALID_FRAC;
    Settings.MinHSubRidgeCoverage = MIN_HSUB_RIDGE_COVERAGE;
    Settings.MinHSubCOIValidFrac = MIN_HSUB_COI_VALID_FRAC;
    Settings.RequireRawPassQC = REQUIRE_RAW_PASSQC;
    Settings.RequireHSubPassQC = REQUIRE_HSUB_PASSQC;
    Settings.URBands = join(UR_BANDS, '; ');
    Settings.HarmonicSensitiveBands = join(HARMONIC_SENSITIVE_BANDS, '; ');

    outXLSX = fullfile(outRoot, 'HSubSupported_PeriodMap.xlsx');
    outMAT  = fullfile(outRoot, 'HSubSupported_PeriodMap.mat');

    validationMap = struct();
    validationMap.Settings = Settings;
    validationMap.LoadSummary = LoadSummary;
    validationMap.AllPeriodCandidates = allCandidates;
    validationMap.Matched_Periods_All = Matched;
    validationMap.CarryForward_Periods = CarryForward;
    validationMap.RawOnly_NotCarriedForward = RawOnly_NotCarriedForward;
    validationMap.HSubOnly_ResidualFeatures = HSubOnly_ResidualFeatures;
    validationMap.FullLadder_Sensitivity = FullLadder_Sensitivity;
    validationMap.Retention_ByBand = Retention_ByBand;
    validationMap.Retention_ByPhotoperiod = Retention_ByPhotoperiod;
    validationMap.Retention_ByPhotoperiodBand = Retention_ByPhotoperiodBand;
    validationMap.QC_Flags = QC_Flags;

    extended_period_gate_write(outXLSX, outMAT, validationMap);

    doFigs = true;
    if isfield(cfg, 'plot') && isfield(cfg.plot, 'generateFigures')
        doFigs = logical(cfg.plot.generateFigures);
    end
    if doFigs
        try
            extended_period_gate_figures(Matched, PrimaryHSubCandidates, figDir, ...
                SAVE_DPI, FIG_EXT, PRIMARY_HSUB_MODE);
        catch ME
            warning('extended_period_gate_run:FiguresFailed', ...
                'Figure generation failed: %s', ME.message);
            extended_period_gate_log(LOG, 'Figure generation failed: %s', ME.message);
        end
    else
        extended_period_gate_log(LOG, 'Figure generation skipped (development / generateFigures=false).');
    end

    extended_period_gate_log(LOG, 'Output workbook: %s', outXLSX);
    extended_period_gate_log(LOG, 'Output MAT: %s', outMAT);
    extended_period_gate_log(LOG, '%s finished at %s', SCRIPT_NAME, ...
        string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    fprintf('\nDone.\n');
    fprintf('Output folder:\n  %s\n', outRoot);
    fprintf('Workbook:\n  %s\n', outXLSX);
    fprintf('MAT file:\n  %s\n', outMAT);
end

%% ------------------------------------------------------------------------
function urNames = resolve_ur_band_names_(cfg)
% Prefer string band names (UR_names). cfg.bands.UR may be numeric limits [n x 2].
    if isfield(cfg, 'bands') && isfield(cfg.bands, 'UR_names') && ~isempty(cfg.bands.UR_names)
        urNames = string(cfg.bands.UR_names(:));
        return;
    end
    if isfield(cfg, 'bands') && isfield(cfg.bands, 'UR') && ~isempty(cfg.bands.UR)
        ur = cfg.bands.UR;
        if isnumeric(ur)
            urNames = ["UR_1_3", "UR_3_6", "UR_6_9", "UR_9_12", "UR_12_18"];
        else
            urNames = string(ur(:));
            urNames = urNames(~ismember(urNames, ["CR_20_28"]));
        end
        return;
    end
    urNames = ["UR_1_3", "UR_3_6", "UR_6_9", "UR_9_12", "UR_12_18"];
end

function cfg = fill_gate_cfg_(cfg)
    if ~isfield(cfg, 'hsub'), cfg.hsub = struct(); end
    if ~isfield(cfg.hsub, 'primaryMode'), cfg.hsub.primaryMode = "SEL_P360"; end
    if ~isfield(cfg.hsub, 'secondaryModes'), cfg.hsub.secondaryModes = ["SEL_P60"]; end
    if ~isfield(cfg.hsub, 'sensitivityModes'), cfg.hsub.sensitivityModes = ["FL_P360", "FL_P60"]; end

    if ~isfield(cfg, 'carryForward'), cfg.carryForward = struct(); end
    if ~isfield(cfg.carryForward, 'periodToleranceFrac'), cfg.carryForward.periodToleranceFrac = 0.15; end
    if ~isfield(cfg.carryForward, 'primaryPhase'), cfg.carryForward.primaryPhase = "All"; end
    if ~isfield(cfg.carryForward, 'minRawRidgeCoverage'), cfg.carryForward.minRawRidgeCoverage = 0.50; end
    if ~isfield(cfg.carryForward, 'minRawCOIValidFrac'), cfg.carryForward.minRawCOIValidFrac = 0.50; end
    if ~isfield(cfg.carryForward, 'minHSubRidgeCoverage'), cfg.carryForward.minHSubRidgeCoverage = 0.50; end
    if ~isfield(cfg.carryForward, 'minHSubCOIValidFrac'), cfg.carryForward.minHSubCOIValidFrac = 0.50; end
    if ~isfield(cfg.carryForward, 'requireRawPassQC'), cfg.carryForward.requireRawPassQC = true; end
    if ~isfield(cfg.carryForward, 'requireHSubPassQC'), cfg.carryForward.requireHSubPassQC = true; end

    if ~isfield(cfg, 'bands'), cfg.bands = struct(); end
    if ~isfield(cfg.bands, 'UR_names') || isempty(cfg.bands.UR_names)
        cfg.bands.UR_names = ["UR_1_3", "UR_3_6", "UR_6_9", "UR_9_12", "UR_12_18"];
    end
    if ~isfield(cfg.bands, 'harmonicSensitive')
        cfg.bands.harmonicSensitive = ["UR_12_18"];
    end

    if ~isfield(cfg, 'plot'), cfg.plot = struct(); end
    if ~isfield(cfg.plot, 'generateFigures'), cfg.plot.generateFigures = true; end
    cfg = extended_apply_plot_cfg(cfg);
end
