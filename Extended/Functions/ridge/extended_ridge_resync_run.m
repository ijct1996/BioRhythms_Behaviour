function resync = extended_ridge_resync_run(handoffDirOrParent, cfg)
%EXTENDED_RIDGE_RESYNC_RUN Ultradian ridge phase resync + BH/FDR + LL projected.
%
%   Modular E2 entry for Kent D (v4 FDR + LL projected aftereffect).
%
%   resync = extended_ridge_resync_run()
%   resync = extended_ridge_resync_run(handoffDir)
%   resync = extended_ridge_resync_run(handoffDir, cfg)
%   resync = extended_ridge_resync_run(opts)   % opts.handoffDir / opts.mapPath / opts.cfg / opts.outputSubdir
%
%   handoffDirOrParent
%     AcrossPhotoperiod_Input folder containing WP_TS__*.mat, OR a parent
%     that contains that folder. Empty → interactive uigetdir / map picker.
%
%   cfg
%     From extended_defaults(); ridge.* and ll.* / plot.* / stats.* fields used.
%
%   Outputs (same as Kent D):
%     {handoffDir}/Ultradian_RidgePhase_Resync/
%       Ultradian_RidgePhase_Resync_Output.xlsx|.mat
%       Figures/, Logs/
%
%   See also: extended_bh_fdr, extended_ll_projected_event_definitions,
%             extended_plot_ridge_power_takeaway, run_extended_ridge_resync

    if nargin < 1, handoffDirOrParent = ''; end
    if nargin < 2, cfg = []; end

    % Allow single-struct options: extended_ridge_resync_run(opts)
    mapPathArg = '';
    outputSubdirArg = '';
    handoffDirExplicit = false;
    if isstruct(handoffDirOrParent) && (isfield(handoffDirOrParent, 'handoffDir') || isfield(handoffDirOrParent, 'mapPath') || isfield(handoffDirOrParent, 'cfg') || isfield(handoffDirOrParent, 'outputSubdir'))
        optsIn = handoffDirOrParent;
        if isfield(optsIn, 'cfg') && ~isempty(optsIn.cfg), cfg = optsIn.cfg; end
        if isfield(optsIn, 'mapPath'), mapPathArg = char(string(optsIn.mapPath)); end
        if isfield(optsIn, 'outputSubdir'), outputSubdirArg = char(string(optsIn.outputSubdir)); end
        if isfield(optsIn, 'handoffDir') && ~isempty(optsIn.handoffDir)
            handoffDirOrParent = optsIn.handoffDir;
            handoffDirExplicit = true;
        else
            handoffDirOrParent = '';
        end
    end

    if isempty(cfg)
        if exist('extended_defaults', 'file') == 2
            cfg = extended_defaults();
        else
            cfg = struct();
        end
    end

    SCRIPT_NAME    = 'Ultradian_RidgePhase_Resync_v4_FDR_LLProjected';
    SCRIPT_VERSION = '4.0-E2';
    RUN_TIMESTAMP  = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

    % ---- Settings (prefer extended_defaults; fall back to Kent D locks) ----
    PRIMARY_PHASE = "All";
    if isfield(cfg, 'carryForward') && isfield(cfg.carryForward, 'primaryPhase')
        PRIMARY_PHASE = string(cfg.carryForward.primaryPhase);
    end
    REQUIRE_VALID_FLAG = true;
    INCLUDE_HARMONIC_SENSITIVE_12H = true;

    UR_BANDS = ["UR_1_3", "UR_3_6", "UR_6_9", "UR_9_12", "UR_12_18"];
    if isfield(cfg, 'bands') && isfield(cfg.bands, 'UR_names') && ~isempty(cfg.bands.UR_names)
        UR_BANDS = string(cfg.bands.UR_names(:));
    elseif isfield(cfg, 'bands') && isfield(cfg.bands, 'UR') && ~isnumeric(cfg.bands.UR)
        UR_BANDS = string(cfg.bands.UR(:));
    end

    PERI_WINDOW_H     = 6.0;
    BIN_WIDTH_H       = 0.5;
    SUMMARY_WINDOW_H  = 2.0;
    MIN_POINTS_PER_BIN = 10;
    MIN_POINTS_SUMMARY = 10;
    if isfield(cfg, 'ridge')
        r = cfg.ridge;
        if isfield(r, 'periWindowH'), PERI_WINDOW_H = r.periWindowH; end
        if isfield(r, 'binWidthH'), BIN_WIDTH_H = r.binWidthH; end
        if isfield(r, 'summaryWindowH'), SUMMARY_WINDOW_H = r.summaryWindowH; end
        if isfield(r, 'minPointsPerBin'), MIN_POINTS_PER_BIN = r.minPointsPerBin; end
        if isfield(r, 'minPointsSummary'), MIN_POINTS_SUMMARY = r.minPointsSummary; end
        if isfield(r, 'requireValidFlag'), REQUIRE_VALID_FLAG = logical(r.requireValidFlag); end
        if isfield(r, 'includeHarmonicSensitive12h'), INCLUDE_HARMONIC_SENSITIVE_12H = logical(r.includeHarmonicSensitive12h); end
    end

    ADAPTIVE_TRANSITION_WINDOWS = true;
    TRANSITION_WINDOW_FRACTION  = 0.45;
    EXCLUDE_INITIAL_DAYS_FROM_RESYNC = false;
    MIN_DAY_FOR_STABLE_RESYNC = 0.0;
    if isfield(cfg, 'ridge')
        r = cfg.ridge;
        if isfield(r, 'adaptiveTransitionWindows'), ADAPTIVE_TRANSITION_WINDOWS = logical(r.adaptiveTransitionWindows); end
        if isfield(r, 'transitionWindowFraction'), TRANSITION_WINDOW_FRACTION = r.transitionWindowFraction; end
        if isfield(r, 'excludeInitialDays'), EXCLUDE_INITIAL_DAYS_FROM_RESYNC = logical(r.excludeInitialDays); end
        if isfield(r, 'minDayForStableResync'), MIN_DAY_FOR_STABLE_RESYNC = r.minDayForStableResync; end
    end

    LL_PHOTOPERIOD_VALUE = 24.0;
    DO_LL_PROJECTED_AFTEREFFECT = true;
    LL_PROJECTED_REFERENCE_PHOTOPERIOD_H = 22.0;
    RUN_LL_PROJECTED_DAY_DECAY = true;
    RUN_L22_VS_LL_PROJECTED_COMPARISON = true;
    if isfield(cfg, 'll')
        ll = cfg.ll;
        if isfield(ll, 'photoperiodValue'), LL_PHOTOPERIOD_VALUE = ll.photoperiodValue; end
        if isfield(ll, 'doProjectedAftereffect'), DO_LL_PROJECTED_AFTEREFFECT = logical(ll.doProjectedAftereffect); end
        if isfield(ll, 'projectedReferencePhotoperiodH'), LL_PROJECTED_REFERENCE_PHOTOPERIOD_H = ll.projectedReferencePhotoperiodH; end
        if isfield(ll, 'runProjectedDayDecay'), RUN_LL_PROJECTED_DAY_DECAY = logical(ll.runProjectedDayDecay); end
        if isfield(ll, 'runL22VsProjectedComparison'), RUN_L22_VS_LL_PROJECTED_COMPARISON = logical(ll.runL22VsProjectedComparison); end
    end

    ALPHA_FDR = 0.05;
    N_PERM_STATS = 10000;
    MIN_CANDIDATES_FOR_TEST = 3;
    if isfield(cfg, 'stats')
        st = cfg.stats;
        if isfield(st, 'alphaFdr'), ALPHA_FDR = st.alphaFdr; end
        if isfield(st, 'nPerm'), N_PERM_STATS = st.nPerm; end
        if isfield(st, 'minCandidatesForTest'), MIN_CANDIDATES_FOR_TEST = st.minCandidatesForTest; end
    end

    SAVE_DPI = 600;
    FIG_EXT  = '.jpg';
    MAX_EXCEL_ROWS_PER_SHEET = 1000000;
    if isfield(cfg, 'plot')
        if isfield(cfg.plot, 'saveDpi'), SAVE_DPI = cfg.plot.saveDpi; end
        if isfield(cfg.plot, 'figExt'), FIG_EXT = char(string(cfg.plot.figExt)); end
    end

    %% ----------------------------- SELECT / RESOLVE INPUTS ------------------
    fprintf('\n%s\n', SCRIPT_NAME);

    if strlength(string(mapPathArg)) > 0 && isfile(mapPathArg)
        mapPath = char(mapPathArg);
        mapDir = fileparts(mapPath);
    else
        fprintf('Select HSubSupported_PeriodMap.mat from the Raw-vs-HSub validation step.\n');
        [mapFile, mapDir] = uigetfile({'*.mat', 'MAT files (*.mat)'}, ...
            'Select HSubSupported_PeriodMap.mat');
        if isequal(mapFile, 0)
            fprintf('No validation map selected. Exiting.\n');
            resync = struct();
            return;
        end
        mapPath = fullfile(mapDir, mapFile);
    end

    Smap = load(mapPath);
    if isfield(Smap, 'validationMap')
        validationMap = Smap.validationMap;
    elseif isfield(Smap, 'CarryForward_Periods')
        validationMap = struct();
        validationMap.CarryForward_Periods = Smap.CarryForward_Periods;
    else
        error('The selected MAT file does not contain validationMap or CarryForward_Periods.');
    end

    if ~isfield(validationMap, 'CarryForward_Periods') || isempty(validationMap.CarryForward_Periods)
        error('No CarryForward_Periods table found in validation map.');
    end
    CarryForward = validationMap.CarryForward_Periods;
    CarryForward = standardise_carryforward_table(CarryForward);

    % Resolve handoffDir: arg → inferred from map → interactive
    handoffDir = '';
    if ~isempty(handoffDirOrParent) && ~(isnumeric(handoffDirOrParent) && isequal(handoffDirOrParent, 0))
        cand = char(string(handoffDirOrParent));
        if ~isempty(dir(fullfile(cand, 'WP_TS__*.mat')))
            handoffDir = cand;
        elseif ~isempty(dir(fullfile(cand, 'AcrossPhotoperiod_Input', 'WP_TS__*.mat')))
            handoffDir = fullfile(cand, 'AcrossPhotoperiod_Input');
        elseif isfolder(cand)
            % Parent of map may already be AcrossPhotoperiod_Input
            handoffDir = cand;
        end
    end

    if isempty(handoffDir) || isempty(dir(fullfile(handoffDir, 'WP_TS__*.mat')))
        inferredHandoffDir = fileparts(mapDir);
        if isempty(dir(fullfile(inferredHandoffDir, 'WP_TS__*.mat')))
            inferredHandoffDir = '';
        end

        if handoffDirExplicit
            if strlength(string(inferredHandoffDir)) > 0
                handoffDir = inferredHandoffDir;
            else
                error('extended_ridge_resync_run:HandoffNotFound', ...
                    'Provided handoffDir has no WP_TS__*.mat and could not infer one from map location: %s', mapPath);
            end
        elseif strlength(string(inferredHandoffDir)) > 0
            useInferred = questdlg(sprintf('Use inferred AcrossPhotoperiod_Input folder?\n\n%s', inferredHandoffDir), ...
                'Use inferred handoff folder?', 'Yes', 'Choose another', 'Yes');
            if strcmpi(useInferred, 'Yes')
                handoffDir = inferredHandoffDir;
            else
                handoffDir = uigetdir(pwd, 'Select AcrossPhotoperiod_Input folder containing WP_TS__*.mat');
            end
        else
            handoffDir = uigetdir(pwd, 'Select AcrossPhotoperiod_Input folder containing WP_TS__*.mat');
        end
    end

    if isequal(handoffDir, 0) || isempty(handoffDir)
        fprintf('No handoff folder selected. Exiting.\n');
        resync = struct();
        return;
    end
    handoffDir = char(handoffDir);

    tsFiles = dir(fullfile(handoffDir, 'WP_TS__*.mat'));
    if isempty(tsFiles)
        error('No WP_TS__*.mat files found in: %s', handoffDir);
    end

    outRoot = fullfile(handoffDir, 'Ultradian_RidgePhase_Resync');
    if strlength(string(outputSubdirArg)) > 0
        outRoot = fullfile(outRoot, outputSubdirArg);
    end
    figDir  = fullfile(outRoot, 'Figures');
    logDir  = fullfile(outRoot, 'Logs');
    ensure_dir(outRoot); ensure_dir(figDir); ensure_dir(logDir);

    logPath = fullfile(logDir, sprintf('%s_Log_%s.txt', SCRIPT_NAME, datestr(now,'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath,'w');
    cleanupObj = onCleanup(@() fclose_if_open(LOG)); %#ok<NASGU>

    log_line(LOG, '%s started at %s', SCRIPT_NAME, RUN_TIMESTAMP);
    log_line(LOG, 'Validation map: %s', mapPath);
    log_line(LOG, 'Handoff folder: %s', handoffDir);
    log_line(LOG, 'Found %d WP_TS files.', numel(tsFiles));

%% ----------------------------- LOAD PHASE TABLES ------------------------
allPhase = table();
loadRows = {};
loadHdr = {'TSMat','Loaded','NRows_RidgePhase','Message'};

for i = 1:numel(tsFiles)
    tsPath = fullfile(tsFiles(i).folder, tsFiles(i).name);
    try
        Sts = load(tsPath);
        if isfield(Sts, 'pkgTS') && isfield(Sts.pkgTS, 'tables') && isfield(Sts.pkgTS.tables, 'RidgePhase_Long')
            T = Sts.pkgTS.tables.RidgePhase_Long;
        elseif isfield(Sts, 'RidgePhase_Long')
            T = Sts.RidgePhase_Long;
        else
            error('No pkgTS.tables.RidgePhase_Long table found.');
        end
        T = standardise_phase_table(T, tsPath);
        allPhase = vertcat_compatible(allPhase, T); %#ok<AGROW>
        loadRows(end+1,:) = {tsPath, true, height(T), 'OK'}; %#ok<SAGROW>
        log_line(LOG, 'Loaded %s (%d RidgePhase rows).', tsPath, height(T));
    catch ME
        loadRows(end+1,:) = {tsPath, false, 0, ME.message}; %#ok<SAGROW>
        log_line(LOG, 'FAILED loading %s: %s', tsPath, ME.message);
    end
end
LoadSummary = cell2table(loadRows, 'VariableNames', loadHdr);

if isempty(allPhase) || height(allPhase) == 0
    error('No RidgePhase_Long rows were loaded. See log: %s', logPath);
end

%% ----------------------------- FILTER TO VALIDATED RAW RIDGES -----------
% Keep only carry-forward rows. These are Raw candidates validated against
% Selective-HSub SEL_P360 by Script 3.
C = CarryForward;
if ismember('CarryForward', C.Properties.VariableNames)
    C = C(logical(C.CarryForward), :);
end
if isempty(C) || height(C) == 0
    error('No carry-forward candidates found in validation map.');
end

if ~INCLUDE_HARMONIC_SENSITIVE_12H && ismember('HarmonicSensitive12hFlag', C.Properties.VariableNames)
    C = C(~logical(C.HarmonicSensitive12hFlag), :);
end

% Prepare candidate lookup table using RawCandidateID as the key.
Ckey = make_candidate_lookup(C);

P = allPhase;
P.SourceNorm = upper(string(P.Source));
P.PhaseNorm  = string(P.Phase);
P.BandNorm   = string(P.BandName);

phaseMask = strcmpi(P.SourceNorm, 'RAW') & strcmpi(P.PhaseNorm, PRIMARY_PHASE) & ismember(P.BandNorm, UR_BANDS);
if REQUIRE_VALID_FLAG && ismember('ValidFlag', P.Properties.VariableNames)
    phaseMask = phaseMask & logical(P.ValidFlag);
end
phaseMask = phaseMask & isfinite(P.RidgePhase_rad) & isfinite(P.RidgePeriod_h) & isfinite(P.Time_days);
P = P(phaseMask, :);

if isempty(P) || height(P) == 0
    error('No valid Raw RidgePhase rows found after filtering.');
end

% Inner join by CandidateID. CandidateID generated by Behav_wavelet_v12 is
% stable and includes file/signal/source/mode/phase/band/rank.
P = innerjoin(P, Ckey, 'Keys', 'CandidateID', 'RightVariables', setdiff(Ckey.Properties.VariableNames, {'CandidateID'}, 'stable'));

if isempty(P) || height(P) == 0
    error('No RidgePhase rows matched the carry-forward validation map. Check CandidateID compatibility.');
end

% v3: optionally remove early days after each photoperiod block starts.
% This is intended to separate stable photoperiod organisation from acute
% adjustment immediately after the weekly 2 h photoperiod shift.
if EXCLUDE_INITIAL_DAYS_FROM_RESYNC && ismember('Time_days', P.Properties.VariableNames)
    nBeforeStableFilter = height(P);
    P = P(P.Time_days >= MIN_DAY_FOR_STABLE_RESYNC, :);
    log_line(LOG, 'Stable-day filter applied to resync phase table: kept %d/%d rows with Time_days >= %.3g.', height(P), nBeforeStableFilter, MIN_DAY_FOR_STABLE_RESYNC);
    if isempty(P) || height(P) == 0
        error('Stable-day filter removed all RidgePhase rows. Lower MIN_DAY_FOR_STABLE_RESYNC or disable EXCLUDE_INITIAL_DAYS_FROM_RESYNC.');
    end
end

log_line(LOG, 'Validated Raw RidgePhase rows after join/filtering: %d', height(P));
log_line(LOG, 'Validated candidates represented: %d', numel(unique(string(P.CandidateID))));

%% ----------------------------- BUILD TRANSITION-LONG TABLE --------------
% v3/v4: LL has no true LD/DL switch and is therefore explicitly summarised
% as a no-transition endpoint. v4 adds a secondary projected-phase analysis
% relative to the former transition phases inherited from the preceding L22 schedule.
LL_NoTransitionSummary = make_ll_no_transition_summary(P, LL_PHOTOPERIOD_VALUE);

TransitionPhase_Long = build_transition_phase_long(P, PERI_WINDOW_H, ADAPTIVE_TRANSITION_WINDOWS, TRANSITION_WINDOW_FRACTION, LL_PHOTOPERIOD_VALUE, LOG);

LL_ProjectedPhase_Long = table();
if DO_LL_PROJECTED_AFTEREFFECT
    LL_ProjectedPhase_Long = build_ll_projected_phase_long(P, PERI_WINDOW_H, ADAPTIVE_TRANSITION_WINDOWS, ...
        TRANSITION_WINDOW_FRACTION, LL_PHOTOPERIOD_VALUE, LL_PROJECTED_REFERENCE_PHOTOPERIOD_H, LOG);
end

if isempty(TransitionPhase_Long) || height(TransitionPhase_Long) == 0
    log_line(LOG, 'No real transition-aligned phase rows were generated. This is expected if only LL/no-transition data were supplied.');
else
    log_line(LOG, 'TransitionPhase_Long rows: %d', height(TransitionPhase_Long));
end
if isempty(LL_ProjectedPhase_Long) || height(LL_ProjectedPhase_Long) == 0
    log_line(LOG, 'No LL projected-transition phase rows were generated.');
else
    log_line(LOG, 'LL_ProjectedPhase_Long rows: %d', height(LL_ProjectedPhase_Long));
end

%% ----------------------------- SUMMARIES --------------------------------
if isempty(TransitionPhase_Long) || height(TransitionPhase_Long) == 0
    BinnedCoherence = table();
    PrePostCoherence = table();
    DeltaR_Summary = table();
    CandidatePrePost = table();
    CandidateDeltaR = table();
else
    BinnedCoherence = make_binned_coherence_summary(TransitionPhase_Long, BIN_WIDTH_H, MIN_POINTS_PER_BIN);
    PrePostCoherence = make_prepost_coherence_summary(TransitionPhase_Long, SUMMARY_WINDOW_H, MIN_POINTS_SUMMARY);
    DeltaR_Summary = make_deltaR_summary(PrePostCoherence);
    CandidatePrePost = make_candidate_prepost_summary(TransitionPhase_Long, SUMMARY_WINDOW_H, MIN_POINTS_SUMMARY);
    CandidateDeltaR = make_candidate_deltaR_summary(CandidatePrePost);
end

% v4 LL projected aftereffect summaries. These are kept separate from real
% LD/DL transition inference because no external transition occurs in LL.
if isempty(LL_ProjectedPhase_Long) || height(LL_ProjectedPhase_Long) == 0
    LL_ProjectedBinnedCoherence = table();
    LL_ProjectedPrePostCoherence = table();
    LL_ProjectedDeltaR = table();
    LL_ProjectedCandidatePrePost = table();
    LL_ProjectedCandidateDeltaR = table();
    LL_ProjectedCandidateDayPrePost = table();
    LL_ProjectedCandidateDayDeltaR = table();
else
    LL_ProjectedBinnedCoherence = make_binned_coherence_summary(LL_ProjectedPhase_Long, BIN_WIDTH_H, MIN_POINTS_PER_BIN);
    LL_ProjectedPrePostCoherence = make_prepost_coherence_summary(LL_ProjectedPhase_Long, SUMMARY_WINDOW_H, MIN_POINTS_SUMMARY);
    LL_ProjectedDeltaR = make_deltaR_summary(LL_ProjectedPrePostCoherence);
    LL_ProjectedCandidatePrePost = make_candidate_prepost_summary(LL_ProjectedPhase_Long, SUMMARY_WINDOW_H, MIN_POINTS_SUMMARY);
    LL_ProjectedCandidateDeltaR = make_candidate_deltaR_summary(LL_ProjectedCandidatePrePost);
    LL_ProjectedCandidateDayPrePost = make_candidate_day_prepost_summary(LL_ProjectedPhase_Long, SUMMARY_WINDOW_H, MIN_POINTS_SUMMARY);
    LL_ProjectedCandidateDayDeltaR = make_candidate_day_deltaR_summary(LL_ProjectedCandidateDayPrePost);
end

% -------------------------- INFERENTIAL STATS + BH/FDR -------------------
% Primary endpoint: candidate-level DeltaR = R_post - R_pre.
% Directional primary test: DeltaR > 0 for real LD/DL transitions.
if isempty(CandidateDeltaR) || height(CandidateDeltaR) == 0
    Resync_PrimaryStats = table();
    Resync_PrimaryStats_BH_FDR = table();
    Resync_RealVsPseudoStats = table();
    Resync_RealVsPseudoStats_BH_FDR = table();
    Resync_RidgePeriodStats = table();
    Resync_RidgePeriodStats_BH_FDR = table();
    Resync_RidgePowerStats = table();
    Resync_RidgePowerStats_BH_FDR = table();
else
    Resync_PrimaryStats = make_resync_primary_stats(CandidateDeltaR, MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    Resync_PrimaryStats_BH_FDR = extended_bh_fdr(Resync_PrimaryStats, 'PValue_raw', ALPHA_FDR);

    % Control endpoint: candidate-level real-transition DeltaR compared against
    % pseudo-transition DeltaR within the same validated candidate where possible.
    Resync_RealVsPseudoStats = make_real_vs_pseudo_stats(CandidateDeltaR, MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    Resync_RealVsPseudoStats_BH_FDR = extended_bh_fdr(Resync_RealVsPseudoStats, 'PValue_raw', ALPHA_FDR);

    % Secondary endpoints: whether ridge period or ridge power changes across the
    % pre/post transition window. These are two-sided and treated separately.
    Resync_RidgePeriodStats = make_candidate_metric_prepost_stats(CandidatePrePost, 'MeanRidgePeriod_h', ...
        'Resync_RidgePeriod_PrePost', 'two-sided', MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    Resync_RidgePeriodStats_BH_FDR = extended_bh_fdr(Resync_RidgePeriodStats, 'PValue_raw', ALPHA_FDR);

    Resync_RidgePowerStats = make_candidate_metric_prepost_stats(CandidatePrePost, 'MeanRidgePower_log10', ...
        'Resync_RidgePower_PrePost', 'two-sided', MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    Resync_RidgePowerStats_BH_FDR = extended_bh_fdr(Resync_RidgePowerStats, 'PValue_raw', ALPHA_FDR);
end

% v4 LL projected aftereffect inference. Kept as a separate FDR family from
% true LD/DL transition resynchronisation.
if isempty(LL_ProjectedCandidateDeltaR) || height(LL_ProjectedCandidateDeltaR) == 0
    LL_ProjectedPrimaryStats = table();
    LL_ProjectedPrimaryStats_BH_FDR = table();
    LL_ProjectedRealVsPseudoStats = table();
    LL_ProjectedRealVsPseudoStats_BH_FDR = table();
    LL_ProjectedRidgePeriodStats = table();
    LL_ProjectedRidgePeriodStats_BH_FDR = table();
    LL_ProjectedRidgePowerStats = table();
    LL_ProjectedRidgePowerStats_BH_FDR = table();
    LL_ProjectDayDecayStats = table();
    LL_ProjectDayDecayStats_BH_FDR = table();
    L22_vs_LLProjectedStats = table();
    L22_vs_LLProjectedStats_BH_FDR = table();
else
    LL_ProjectedPrimaryStats = make_projected_primary_stats(LL_ProjectedCandidateDeltaR, MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    LL_ProjectedPrimaryStats_BH_FDR = extended_bh_fdr(LL_ProjectedPrimaryStats, 'PValue_raw', ALPHA_FDR);

    LL_ProjectedRealVsPseudoStats = make_projected_real_vs_pseudo_stats(LL_ProjectedCandidateDeltaR, MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    LL_ProjectedRealVsPseudoStats_BH_FDR = extended_bh_fdr(LL_ProjectedRealVsPseudoStats, 'PValue_raw', ALPHA_FDR);

    LL_ProjectedRidgePeriodStats = make_candidate_metric_prepost_stats(LL_ProjectedCandidatePrePost, 'MeanRidgePeriod_h', ...
        'LL_Projected_RidgePeriod_PrePost', 'two-sided', MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    LL_ProjectedRidgePeriodStats_BH_FDR = extended_bh_fdr(LL_ProjectedRidgePeriodStats, 'PValue_raw', ALPHA_FDR);

    LL_ProjectedRidgePowerStats = make_candidate_metric_prepost_stats(LL_ProjectedCandidatePrePost, 'MeanRidgePower_log10', ...
        'LL_Projected_RidgePower_PrePost', 'two-sided', MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
    LL_ProjectedRidgePowerStats_BH_FDR = extended_bh_fdr(LL_ProjectedRidgePowerStats, 'PValue_raw', ALPHA_FDR);

    if RUN_LL_PROJECTED_DAY_DECAY
        LL_ProjectDayDecayStats = make_ll_projected_day_decay_stats(LL_ProjectedCandidateDayDeltaR, MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
        LL_ProjectDayDecayStats_BH_FDR = extended_bh_fdr(LL_ProjectDayDecayStats, 'PValue_raw', ALPHA_FDR);
    else
        LL_ProjectDayDecayStats = table();
        LL_ProjectDayDecayStats_BH_FDR = table();
    end

    if RUN_L22_VS_LL_PROJECTED_COMPARISON
        L22_vs_LLProjectedStats = make_l22_vs_ll_projected_stats(CandidateDeltaR, LL_ProjectedCandidateDeltaR, ...
            LL_PROJECTED_REFERENCE_PHOTOPERIOD_H, MIN_CANDIDATES_FOR_TEST, N_PERM_STATS);
        L22_vs_LLProjectedStats_BH_FDR = extended_bh_fdr(L22_vs_LLProjectedStats, 'PValue_raw', ALPHA_FDR);
    else
        L22_vs_LLProjectedStats = table();
        L22_vs_LLProjectedStats_BH_FDR = table();
    end
end

Settings = make_settings_table(SCRIPT_NAME, SCRIPT_VERSION, RUN_TIMESTAMP, mapPath, handoffDir, ...
    PRIMARY_PHASE, REQUIRE_VALID_FLAG, INCLUDE_HARMONIC_SENSITIVE_12H, UR_BANDS, ...
    PERI_WINDOW_H, BIN_WIDTH_H, SUMMARY_WINDOW_H, MIN_POINTS_PER_BIN, MIN_POINTS_SUMMARY, ...
    ALPHA_FDR, N_PERM_STATS, MIN_CANDIDATES_FOR_TEST, ...
    ADAPTIVE_TRANSITION_WINDOWS, TRANSITION_WINDOW_FRACTION, EXCLUDE_INITIAL_DAYS_FROM_RESYNC, MIN_DAY_FOR_STABLE_RESYNC, LL_PHOTOPERIOD_VALUE);
Settings.DoLLProjectedAftereffect = logical(DO_LL_PROJECTED_AFTEREFFECT);
Settings.LLProjectedReferencePhotoperiod_h = LL_PROJECTED_REFERENCE_PHOTOPERIOD_H;
Settings.RunLLProjectedDayDecay = logical(RUN_LL_PROJECTED_DAY_DECAY);
Settings.RunL22VsLLProjectedComparison = logical(RUN_L22_VS_LL_PROJECTED_COMPARISON);
Settings.LLProjectedHandling = "Secondary projected-phase aftereffect analysis only; LL is not treated as having real LD/DL transitions";

%% ----------------------------- WRITE OUTPUTS ----------------------------
outXLSX = fullfile(outRoot, 'Ultradian_RidgePhase_Resync_Output.xlsx');
outMAT  = fullfile(outRoot, 'Ultradian_RidgePhase_Resync_Output.mat');
delete_if_exists(outXLSX);

writecell({'Ridge-following phase resynchronisation analysis.', ...
           'Input phase source: Raw RidgePhase_Long rows only.', ...
           'Candidate inclusion: CarryForward_Periods from Raw-vs-Selective-HSub validation.', ...
           'DeltaR = R_post - R_pre, where R is circular resultant length.', ...
           'Controls: MidLight and MidDark pseudo-transitions.', ...
           'Inferential p-values are candidate-level sign-flip permutation tests.', ...
           'BH/FDR correction is applied separately within predefined analysis families.', ...
           'LL/no-transition data are not analysed as real LD/DL transitions; they are reported separately in LL_NoTransitionSummary.', ...
           'Secondary LL analysis: LL ridge phase is aligned to projected/former transition phases inherited from the preceding L22 schedule.', ...
           'LL projected event labels are ProjectedDL_LL, ProjectedLD_LL, ProjectedMidLight_LL and ProjectedMidDark_LL.', ...
           'For short-dark photoperiods, transition windows are adaptively capped to avoid overlap across neighbouring transitions.'}', outXLSX, 'Sheet', 'README');

safe_writetable_xlsx(Settings, outXLSX, 'Settings', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LoadSummary, outXLSX, 'LoadSummary', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Ckey, outXLSX, 'ValidatedCandidates_Used', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_NoTransitionSummary, outXLSX, 'LL_NoTransitionSummary', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedPhase_Long, outXLSX, 'LL_ProjectedPhase_Long', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedBinnedCoherence, outXLSX, 'LLProj_BinnedCoherence', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedPrePostCoherence, outXLSX, 'LLProj_PrePostCoherence', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedDeltaR, outXLSX, 'LLProj_DeltaR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedCandidatePrePost, outXLSX, 'LLProj_CandidatePrePost', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedCandidateDeltaR, outXLSX, 'LLProj_CandidateDeltaR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedCandidateDayPrePost, outXLSX, 'LLProj_CandDayPrePost', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedCandidateDayDeltaR, outXLSX, 'LLProj_CandDayDeltaR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedPrimaryStats, outXLSX, 'LLProj_PrimaryStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedPrimaryStats_BH_FDR, outXLSX, 'LLProj_PrimaryStats_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedRealVsPseudoStats, outXLSX, 'LLProj_RealVsPseudoStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedRealVsPseudoStats_BH_FDR, outXLSX, 'LLProj_RealVsPseudo_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedRidgePeriodStats, outXLSX, 'LLProj_RidgePeriodStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedRidgePeriodStats_BH_FDR, outXLSX, 'LLProj_RidgePeriod_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedRidgePowerStats, outXLSX, 'LLProj_RidgePowerStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectedRidgePowerStats_BH_FDR, outXLSX, 'LLProj_RidgePower_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectDayDecayStats, outXLSX, 'LLProj_DayDecayStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(LL_ProjectDayDecayStats_BH_FDR, outXLSX, 'LLProj_DayDecay_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(L22_vs_LLProjectedStats, outXLSX, 'L22_vs_LLProjectedStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(L22_vs_LLProjectedStats_BH_FDR, outXLSX, 'L22_vs_LLProj_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(TransitionPhase_Long, outXLSX, 'TransitionPhase_Long', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(BinnedCoherence, outXLSX, 'BinnedCoherence', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(PrePostCoherence, outXLSX, 'PrePostCoherence', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(DeltaR_Summary, outXLSX, 'DeltaR_Summary', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(CandidatePrePost, outXLSX, 'CandidatePrePost', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(CandidateDeltaR, outXLSX, 'CandidateDeltaR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_PrimaryStats, outXLSX, 'Resync_PrimaryStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_PrimaryStats_BH_FDR, outXLSX, 'Resync_PrimaryStats_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_RealVsPseudoStats, outXLSX, 'Resync_RealVsPseudoStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_RealVsPseudoStats_BH_FDR, outXLSX, 'RealVsPseudoStats_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_RidgePeriodStats, outXLSX, 'RidgePeriodStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_RidgePeriodStats_BH_FDR, outXLSX, 'RidgePeriodStats_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_RidgePowerStats, outXLSX, 'RidgePowerStats', outRoot, MAX_EXCEL_ROWS_PER_SHEET);
safe_writetable_xlsx(Resync_RidgePowerStats_BH_FDR, outXLSX, 'RidgePowerStats_BH_FDR', outRoot, MAX_EXCEL_ROWS_PER_SHEET);

resync = struct();
resync.Settings = Settings;
resync.LoadSummary = LoadSummary;
resync.ValidatedCandidates_Used = Ckey;
resync.LL_NoTransitionSummary = LL_NoTransitionSummary;
resync.LL_ProjectedPhase_Long = LL_ProjectedPhase_Long;
resync.LL_ProjectedBinnedCoherence = LL_ProjectedBinnedCoherence;
resync.LL_ProjectedPrePostCoherence = LL_ProjectedPrePostCoherence;
resync.LL_ProjectedDeltaR = LL_ProjectedDeltaR;
resync.LL_ProjectedCandidatePrePost = LL_ProjectedCandidatePrePost;
resync.LL_ProjectedCandidateDeltaR = LL_ProjectedCandidateDeltaR;
resync.LL_ProjectedCandidateDayPrePost = LL_ProjectedCandidateDayPrePost;
resync.LL_ProjectedCandidateDayDeltaR = LL_ProjectedCandidateDayDeltaR;
resync.LL_ProjectedPrimaryStats = LL_ProjectedPrimaryStats;
resync.LL_ProjectedPrimaryStats_BH_FDR = LL_ProjectedPrimaryStats_BH_FDR;
resync.LL_ProjectedRealVsPseudoStats = LL_ProjectedRealVsPseudoStats;
resync.LL_ProjectedRealVsPseudoStats_BH_FDR = LL_ProjectedRealVsPseudoStats_BH_FDR;
resync.LL_ProjectedRidgePeriodStats = LL_ProjectedRidgePeriodStats;
resync.LL_ProjectedRidgePeriodStats_BH_FDR = LL_ProjectedRidgePeriodStats_BH_FDR;
resync.LL_ProjectedRidgePowerStats = LL_ProjectedRidgePowerStats;
resync.LL_ProjectedRidgePowerStats_BH_FDR = LL_ProjectedRidgePowerStats_BH_FDR;
resync.LL_ProjectDayDecayStats = LL_ProjectDayDecayStats;
resync.LL_ProjectDayDecayStats_BH_FDR = LL_ProjectDayDecayStats_BH_FDR;
resync.L22_vs_LLProjectedStats = L22_vs_LLProjectedStats;
resync.L22_vs_LLProjectedStats_BH_FDR = L22_vs_LLProjectedStats_BH_FDR;
resync.TransitionPhase_Long = TransitionPhase_Long;
resync.BinnedCoherence = BinnedCoherence;
resync.PrePostCoherence = PrePostCoherence;
resync.DeltaR_Summary = DeltaR_Summary;
resync.CandidatePrePost = CandidatePrePost;
resync.CandidateDeltaR = CandidateDeltaR;
resync.Resync_PrimaryStats = Resync_PrimaryStats;
resync.Resync_PrimaryStats_BH_FDR = Resync_PrimaryStats_BH_FDR;
resync.Resync_RealVsPseudoStats = Resync_RealVsPseudoStats;
resync.Resync_RealVsPseudoStats_BH_FDR = Resync_RealVsPseudoStats_BH_FDR;
resync.Resync_RidgePeriodStats = Resync_RidgePeriodStats;
resync.Resync_RidgePeriodStats_BH_FDR = Resync_RidgePeriodStats_BH_FDR;
resync.Resync_RidgePowerStats = Resync_RidgePowerStats;
resync.Resync_RidgePowerStats_BH_FDR = Resync_RidgePowerStats_BH_FDR;
save(outMAT, 'resync', '-v7.3');

%% ----------------------------- FIGURES ----------------------------------
doFigs = true;
if isfield(cfg, 'plot') && isfield(cfg.plot, 'generateFigures')
    doFigs = logical(cfg.plot.generateFigures);
end
if doFigs
    try
        make_resync_figures(BinnedCoherence, DeltaR_Summary, CandidateDeltaR, figDir, SAVE_DPI, FIG_EXT);
        make_ll_projected_figures(LL_ProjectedBinnedCoherence, LL_ProjectedDeltaR, LL_ProjectedCandidateDeltaR, figDir, SAVE_DPI, FIG_EXT);
    catch ME
        warning('Figure generation failed: %s', ME.message);
        log_line(LOG, 'Figure generation failed: %s', ME.message);
    end
else
    log_line(LOG, 'Resync figures skipped (generateFigures=false).');
end

log_line(LOG, 'Output workbook: %s', outXLSX);
log_line(LOG, 'Output MAT: %s', outMAT);
log_line(LOG, '%s finished at %s', SCRIPT_NAME, string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));

fprintf('\nDone.\n');
fprintf('Output folder:\n  %s\n', outRoot);

    resync.outRoot = outRoot;
    resync.outXLSX = outXLSX;
    resync.outMAT = outMAT;
    resync.handoffDir = handoffDir;
    resync.mapPath = mapPath;
end

function E = projected_ll_event_definitions(refPhotoperiodH)
    % Compatibility shim → Extended/Functions/ridge/
    E = extended_ll_projected_event_definitions(refPhotoperiodH);
end

function T = add_bh_fdr(T, pVar, alpha)
    % Compatibility shim → Extended/Functions/stats/
    T = extended_bh_fdr(T, pVar, alpha);
end
function T = standardise_carryforward_table(T)
    if ~istable(T), error('CarryForward_Periods is not a table.'); end
    strVars = {'File','SignalID','ConditionParsed','Phase','BandName','RawCandidateID', ...
               'PrimaryHSubMode','PrimaryHSubCandidateID','MatchStatus','QCFlag', ...
               'FullLadderSensitivityStatus','FinalValidationClass','SourcePackage'};
    for i = 1:numel(strVars)
        v = strVars{i};
        if ismember(v, T.Properties.VariableNames)
            T.(v) = string(T.(v));
        end
    end
    numVars = {'Photoperiod_h','RawPeriod_h','RawIQR_h','RawMeanBandPower_log10','RawMeanRidgePower_log10', ...
               'RawRidgeCoverageFrac','RawCOIValidFrac','PrimaryHSubPeriod_h','PeriodDiffPercent'};
    for i = 1:numel(numVars)
        v = numVars{i};
        if ismember(v, T.Properties.VariableNames)
            T.(v) = double(T.(v));
        end
    end
    if ismember('CarryForward', T.Properties.VariableNames)
        T.CarryForward = to_logical(T.CarryForward);
    end
    if ismember('HarmonicSensitive12hFlag', T.Properties.VariableNames)
        T.HarmonicSensitive12hFlag = to_logical(T.HarmonicSensitive12hFlag);
    else
        T.HarmonicSensitive12hFlag = false(height(T),1);
    end
end

function T = standardise_phase_table(T, sourceMat)
    if ~istable(T), error('RidgePhase_Long object is not a table.'); end
    required = {'File','SignalID','ConditionParsed','Source','HSubResidualMode','Photoperiod_h','BandName','CandidateID', ...
                'Time_days','ZT_hours','LightStateValue','Phase','RidgePeriod_h','RidgePower_log10','RidgePhase_rad','ValidFlag'};
    for i = 1:numel(required)
        if ~ismember(required{i}, T.Properties.VariableNames)
            error('RidgePhase_Long is missing required column: %s', required{i});
        end
    end
    strVars = {'File','SignalID','ConditionParsed','Source','HSubResidualMode','BandName','CandidateID','LightStateValue','Phase'};
    for i = 1:numel(strVars)
        T.(strVars{i}) = string(T.(strVars{i}));
    end
    numVars = {'Photoperiod_h','Time_days','ZT_hours','RidgePeriod_h','RidgePower_log10','RidgePhase_rad'};
    for i = 1:numel(numVars)
        T.(numVars{i}) = double(T.(numVars{i}));
    end
    T.ValidFlag = to_logical(T.ValidFlag);
    T.SourcePackage = repmat(string(sourceMat), height(T), 1);
end

function Ckey = make_candidate_lookup(C)
    needed = {'RawCandidateID','File','SignalID','ConditionParsed','Photoperiod_h','Phase','BandName','RawPeriod_h'};
    for i = 1:numel(needed)
        if ~ismember(needed{i}, C.Properties.VariableNames)
            error('CarryForward_Periods is missing required column: %s', needed{i});
        end
    end
    Ckey = table();
    Ckey.CandidateID = string(C.RawCandidateID);
    Ckey.File_validated = string(C.File);
    Ckey.SignalID_validated = string(C.SignalID);
    Ckey.ConditionParsed_validated = string(C.ConditionParsed);
    Ckey.Photoperiod_h_validated = double(C.Photoperiod_h);
    Ckey.Phase_validated = string(C.Phase);
    Ckey.BandName_validated = string(C.BandName);
    Ckey.RawPeriod_h_validated = double(C.RawPeriod_h);

    optionalMap = {'RawIQR_h','RawIQR_h_validated'; ...
                   'PeriodDiffPercent','RawVsHSubPeriodDiffPercent'; ...
                   'PrimaryHSubPeriod_h','PrimaryHSubPeriod_h_validated'; ...
                   'RawRidgeCoverageFrac','RawRidgeCoverageFrac_validated'; ...
                   'RawCOIValidFrac','RawCOIValidFrac_validated'; ...
                   'FullLadderSensitivityStatus','FullLadderSensitivityStatus'; ...
                   'FinalValidationClass','FinalValidationClass'; ...
                   'HarmonicSensitive12hFlag','HarmonicSensitive12hFlag'};
    for i = 1:size(optionalMap,1)
        src = optionalMap{i,1}; dst = optionalMap{i,2};
        if ismember(src, C.Properties.VariableNames)
            Ckey.(dst) = C.(src);
        end
    end

    % Ensure one row per CandidateID.
    [~, ia] = unique(Ckey.CandidateID, 'stable');
    Ckey = Ckey(ia,:);
end

function TL = build_transition_phase_long(P, periWindowH, adaptiveWindows, windowFraction, llPhotoValue, LOG)
    vars = {'File','SignalID','ConditionParsed','Photoperiod_h','BandName','CandidateID', ...
            'RawPeriod_h_validated','PrimaryHSubPeriod_h_validated','RawVsHSubPeriodDiffPercent', ...
            'FullLadderSensitivityStatus','FinalValidationClass','HarmonicSensitive12hFlag', ...
            'TransitionType','TransitionZT_h','TransitionDay','TransitionAbsTime_h', ...
            'EffectivePreWindow_h','EffectivePostWindow_h', ...
            'RelativeTime_h','RelativeCycles','Time_days','ZT_hours','LightStateValue', ...
            'RidgePhase_rad','RidgePeriod_h','RidgePower_log10','ValidFlag'};
    rows = {};

    candIDs = unique(string(P.CandidateID), 'stable');
    log_line(LOG, 'Building transition windows for %d validated candidates.', numel(candIDs));

    for c = 1:numel(candIDs)
        Pc = P(string(P.CandidateID) == candIDs(c), :);
        if isempty(Pc), continue; end

        Pc = sortrows(Pc, 'Time_days');
        photo = robust_scalar_mode(Pc.Photoperiod_h);
        if ~isfinite(photo) || photo <= 0
            photo = robust_scalar_mode(Pc.Photoperiod_h_validated);
        end
        if ~isfinite(photo) || photo <= 0
            log_line(LOG, 'Candidate %s skipped: invalid photoperiod.', candIDs(c));
            continue;
        end
        if photo >= llPhotoValue || photo >= 24
            log_line(LOG, 'Candidate %s treated as LL/no-transition and excluded from LD/DL transition windows.', candIDs(c));
            continue;
        end

        rawPeriod = robust_scalar_mode(Pc.RawPeriod_h_validated);
        if ~isfinite(rawPeriod) || rawPeriod <= 0
            rawPeriod = median(Pc.RidgePeriod_h, 'omitnan');
        end
        if ~isfinite(rawPeriod) || rawPeriod <= 0
            continue;
        end

        tH = Pc.Time_days * 24;
        tMin = min(tH, [], 'omitnan');
        tMax = max(tH, [], 'omitnan');
        if ~isfinite(tMin) || ~isfinite(tMax), continue; end

        events = transition_definitions(photo);
        dayStart = floor(tMin/24) - 1;
        dayEnd   = ceil(tMax/24) + 1;

        for e = 1:height(events)
            evType = string(events.TransitionType(e));
            evZT   = events.TransitionZT_h(e);
            [preWin, postWin] = transition_window_limits(photo, evType, periWindowH, adaptiveWindows, windowFraction);
            if ~isfinite(preWin) || ~isfinite(postWin) || preWin <= 0 || postWin <= 0
                continue;
            end

            for d = dayStart:dayEnd
                evAbs = d*24 + evZT;
                if evAbs < (tMin - preWin) || evAbs > (tMax + postWin)
                    continue;
                end
                rel = tH - evAbs;
                idx = rel >= -preWin & rel <= postWin & isfinite(Pc.RidgePhase_rad) & logical(Pc.ValidFlag);
                if ~any(idx), continue; end

                Pi = Pc(idx,:);
                rel_i = rel(idx);
                relCycles_i = rel_i ./ rawPeriod;
                nAdd = height(Pi);
                for r = 1:nAdd
                    rows(end+1,:) = {string(Pi.File(r)), string(Pi.SignalID(r)), string(Pi.ConditionParsed(r)), double(Pi.Photoperiod_h(r)), ...
                        string(Pi.BandName(r)), string(Pi.CandidateID(r)), double(Pi.RawPeriod_h_validated(r)), get_optional_num(Pi,'PrimaryHSubPeriod_h_validated',r), ...
                        get_optional_num(Pi,'RawVsHSubPeriodDiffPercent',r), get_optional_str(Pi,'FullLadderSensitivityStatus',r), ...
                        get_optional_str(Pi,'FinalValidationClass',r), get_optional_logical(Pi,'HarmonicSensitive12hFlag',r), ...
                        evType, evZT, d, evAbs, preWin, postWin, rel_i(r), relCycles_i(r), double(Pi.Time_days(r)), double(Pi.ZT_hours(r)), ...
                        string(Pi.LightStateValue(r)), double(Pi.RidgePhase_rad(r)), double(Pi.RidgePeriod_h(r)), double(Pi.RidgePower_log10(r)), logical(Pi.ValidFlag(r))}; %#ok<AGROW>
                end
            end
        end
    end

    if isempty(rows)
        TL = cell2table(cell(0,numel(vars)), 'VariableNames', vars);
    else
        TL = cell2table(rows, 'VariableNames', vars);
    end
end


function TL = build_ll_projected_phase_long(P, periWindowH, adaptiveWindows, windowFraction, llPhotoValue, refPhotoperiodH, LOG)
    % v4: LL projected aftereffect analysis.
    % LL is aligned to former/projected transition phases inherited from the
    % preceding LD schedule. This does not define real light transitions in LL.
    vars = {'File','SignalID','ConditionParsed','Photoperiod_h','BandName','CandidateID', ...
            'RawPeriod_h_validated','PrimaryHSubPeriod_h_validated','RawVsHSubPeriodDiffPercent', ...
            'FullLadderSensitivityStatus','FinalValidationClass','HarmonicSensitive12hFlag', ...
            'TransitionType','TransitionZT_h','TransitionDay','TransitionAbsTime_h', ...
            'EffectivePreWindow_h','EffectivePostWindow_h', ...
            'RelativeTime_h','RelativeCycles','Time_days','ZT_hours','LightStateValue', ...
            'RidgePhase_rad','RidgePeriod_h','RidgePower_log10','ValidFlag'};
    rows = {};

    if isempty(P) || height(P)==0
        TL = cell2table(cell(0,numel(vars)), 'VariableNames', vars); return;
    end

    photo = P.Photoperiod_h;
    if ismember('Photoperiod_h_validated', P.Properties.VariableNames)
        idxBad = ~isfinite(photo);
        photo(idxBad) = P.Photoperiod_h_validated(idxBad);
    end
    PLL = P(isfinite(photo) & (photo >= llPhotoValue | photo >= 24), :);
    if isempty(PLL)
        TL = cell2table(cell(0,numel(vars)), 'VariableNames', vars); return;
    end

    candIDs = unique(string(PLL.CandidateID), 'stable');
    log_line(LOG, 'Building projected LL windows for %d validated LL candidates using former L%.3g schedule.', numel(candIDs), refPhotoperiodH);

    events = extended_ll_projected_event_definitions(refPhotoperiodH);

    for c = 1:numel(candIDs)
        Pc = PLL(string(PLL.CandidateID) == candIDs(c), :);
        if isempty(Pc), continue; end
        Pc = sortrows(Pc, 'Time_days');

        rawPeriod = robust_scalar_mode(Pc.RawPeriod_h_validated);
        if ~isfinite(rawPeriod) || rawPeriod <= 0
            rawPeriod = median(Pc.RidgePeriod_h, 'omitnan');
        end
        if ~isfinite(rawPeriod) || rawPeriod <= 0
            continue;
        end

        tH = Pc.Time_days * 24;
        tMin = min(tH, [], 'omitnan');
        tMax = max(tH, [], 'omitnan');
        if ~isfinite(tMin) || ~isfinite(tMax), continue; end

        dayStart = floor(tMin/24) - 1;
        dayEnd   = ceil(tMax/24) + 1;

        for e = 1:height(events)
            evType = string(events.TransitionType(e));
            evZT   = events.TransitionZT_h(e);
            baseType = string(events.BaseTransitionType(e));
            [preWin, postWin] = projected_transition_window_limits(refPhotoperiodH, baseType, periWindowH, adaptiveWindows, windowFraction);
            if ~isfinite(preWin) || ~isfinite(postWin) || preWin <= 0 || postWin <= 0
                continue;
            end

            for d = dayStart:dayEnd
                evAbs = d*24 + evZT;
                if evAbs < (tMin - preWin) || evAbs > (tMax + postWin)
                    continue;
                end
                rel = tH - evAbs;
                idx = rel >= -preWin & rel <= postWin & isfinite(Pc.RidgePhase_rad) & logical(Pc.ValidFlag);
                if ~any(idx), continue; end

                Pi = Pc(idx,:);
                rel_i = rel(idx);
                relCycles_i = rel_i ./ rawPeriod;
                nAdd = height(Pi);
                for r = 1:nAdd
                    rows(end+1,:) = {string(Pi.File(r)), string(Pi.SignalID(r)), string(Pi.ConditionParsed(r)), double(Pi.Photoperiod_h(r)), ...
                        string(Pi.BandName(r)), string(Pi.CandidateID(r)), double(Pi.RawPeriod_h_validated(r)), get_optional_num(Pi,'PrimaryHSubPeriod_h_validated',r), ...
                        get_optional_num(Pi,'RawVsHSubPeriodDiffPercent',r), get_optional_str(Pi,'FullLadderSensitivityStatus',r), ...
                        get_optional_str(Pi,'FinalValidationClass',r), get_optional_logical(Pi,'HarmonicSensitive12hFlag',r), ...
                        evType, evZT, d, evAbs, preWin, postWin, rel_i(r), relCycles_i(r), double(Pi.Time_days(r)), double(Pi.ZT_hours(r)), ...
                        string(Pi.LightStateValue(r)), double(Pi.RidgePhase_rad(r)), double(Pi.RidgePeriod_h(r)), double(Pi.RidgePower_log10(r)), logical(Pi.ValidFlag(r))}; %#ok<AGROW>
                end
            end
        end
    end

    if isempty(rows)
        TL = cell2table(cell(0,numel(vars)), 'VariableNames', vars);
    else
        TL = cell2table(rows, 'VariableNames', vars);
    end
end

function [preWin, postWin] = projected_transition_window_limits(refPhotoperiodH, baseTransitionType, periWindowH, adaptiveWindows, windowFraction)
    [preWin, postWin] = transition_window_limits(refPhotoperiodH, baseTransitionType, periWindowH, adaptiveWindows, windowFraction);
end

function [preWin, postWin] = transition_window_limits(photoperiod_h, transitionType, periWindowH, adaptiveWindows, windowFraction)
    % v3: protect compressed-dark photoperiods (e.g., L20/L22) by preventing
    % peri-transition windows from extending deep into the neighbouring transition.
    photo = double(photoperiod_h);
    dark_h = 24 - photo;
    preWin = periWindowH;
    postWin = periWindowH;

    if ~adaptiveWindows
        return;
    end

    transitionType = string(transitionType);
    switch transitionType
        case "DL"       % lights-on: previous interval is dark, next interval is light
            preCap  = windowFraction * dark_h;
            postCap = windowFraction * photo;
        case "LD"       % lights-off: previous interval is light, next interval is dark
            preCap  = windowFraction * photo;
            postCap = windowFraction * dark_h;
        case "MidLight"
            preCap  = windowFraction * (photo/2);
            postCap = windowFraction * (photo/2);
        case "MidDark"
            preCap  = windowFraction * (dark_h/2);
            postCap = windowFraction * (dark_h/2);
        otherwise
            preCap = periWindowH;
            postCap = periWindowH;
    end

    preWin  = min(periWindowH, max(0, preCap));
    postWin = min(periWindowH, max(0, postCap));
end

function LL = make_ll_no_transition_summary(P, llPhotoValue)
    hdr = {'File','SignalID','ConditionParsed','Photoperiod_h','BandName','CandidateID', ...
           'ReasonExcludedFromTransitionAnalysis','N_PhaseRows','TimeStart_day','TimeEnd_day', ...
           'MeanRidgePeriod_h','SDRidgePeriod_h','MeanRidgePower_log10','SDRidgePower_log10', ...
           'R_AllRecording','MeanPhase_AllRecording_rad','CircularSD_AllRecording_rad'};
    rows = {};
    if isempty(P) || height(P) == 0
        LL = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr);
        return;
    end

    photo = P.Photoperiod_h;
    if ismember('Photoperiod_h_validated', P.Properties.VariableNames)
        idxBad = ~isfinite(photo);
        photo(idxBad) = P.Photoperiod_h_validated(idxBad);
    end
    isLL = isfinite(photo) & (photo >= llPhotoValue | photo >= 24);
    if ~any(isLL)
        LL = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr);
        return;
    end

    PL = P(isLL,:);
    groups = unique(PL(:, {'File','SignalID','ConditionParsed','Photoperiod_h','BandName','CandidateID'}), 'rows', 'stable');
    for g = 1:height(groups)
        idx = string(PL.File)==string(groups.File(g)) & string(PL.SignalID)==string(groups.SignalID(g)) & ...
              string(PL.ConditionParsed)==string(groups.ConditionParsed(g)) & PL.Photoperiod_h==groups.Photoperiod_h(g) & ...
              string(PL.BandName)==string(groups.BandName(g)) & string(PL.CandidateID)==string(groups.CandidateID(g));
        G = PL(idx,:);
        theta = G.RidgePhase_rad(isfinite(G.RidgePhase_rad));
        [R, mu, csd] = circ_summary(theta);
        rows(end+1,:) = {string(groups.File(g)), string(groups.SignalID(g)), string(groups.ConditionParsed(g)), double(groups.Photoperiod_h(g)), ...
            string(groups.BandName(g)), string(groups.CandidateID(g)), ...
            "LL/no true LD or DL transition; excluded from transition-resynchronisation inference", ...
            height(G), min(G.Time_days,[],'omitnan'), max(G.Time_days,[],'omitnan'), ...
            mean(G.RidgePeriod_h,'omitnan'), std(G.RidgePeriod_h,'omitnan'), ...
            mean(G.RidgePower_log10,'omitnan'), std(G.RidgePower_log10,'omitnan'), R, mu, csd}; %#ok<AGROW>
    end
    LL = rows_to_table(rows, hdr);
end

function E = transition_definitions(photoperiod_h)
    dark_h = 24 - photoperiod_h;
    types = ["DL"; "LD"; "MidLight"; "MidDark"];
    zts = [0; photoperiod_h; photoperiod_h/2; mod(photoperiod_h + dark_h/2, 24)];
    E = table(types, zts, 'VariableNames', {'TransitionType','TransitionZT_h'});
end

function S = make_binned_coherence_summary(TL, binWidthH, minN)
    edges = -ceil(max(abs(TL.RelativeTime_h))) : binWidthH : ceil(max(abs(TL.RelativeTime_h)));
    if numel(edges) < 2
        edges = -6:binWidthH:6;
    end
    centers = edges(1:end-1) + diff(edges)/2;
    binIdx = discretize(TL.RelativeTime_h, edges);
    TL.RelBin = binIdx;

    groups = unique(TL(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    rows = {};
    hdr = {'Photoperiod_h','BandName','TransitionType','RelBin','RelBinCenter_h','RelBinStart_h','RelBinEnd_h', ...
           'N_PhaseObs','N_Mice','N_Candidates','N_Transitions','R','MeanPhase_rad','CircularSD_rad', ...
           'MeanRidgePeriod_h','MeanRidgePower_log10','PassN'};

    for g = 1:height(groups)
        gm = TL.Photoperiod_h == groups.Photoperiod_h(g) & ...
             string(TL.BandName) == string(groups.BandName(g)) & ...
             string(TL.TransitionType) == string(groups.TransitionType(g));
        for b = 1:numel(centers)
            idx = gm & TL.RelBin == b;
            theta = TL.RidgePhase_rad(idx);
            theta = theta(isfinite(theta));
            n = numel(theta);
            [R, mu, csd] = circ_summary(theta);
            passN = n >= minN;
            rows(end+1,:) = {groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), b, centers(b), edges(b), edges(b+1), ...
                n, count_unique(TL.SignalID(idx)), count_unique(TL.CandidateID(idx)), count_unique(string(TL.CandidateID(idx)) + "|" + string(TL.TransitionType(idx)) + "|" + string(TL.TransitionDay(idx))), ...
                R, mu, csd, mean(TL.RidgePeriod_h(idx),'omitnan'), mean(TL.RidgePower_log10(idx),'omitnan'), passN}; %#ok<AGROW>
        end
    end

    S = rows_to_table(rows, hdr);
end

function S = make_prepost_coherence_summary(TL, summaryWindowH, minN)
    groups = unique(TL(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    winNames = ["Pre"; "Post"];
    rows = {};
    hdr = {'Photoperiod_h','BandName','TransitionType','Window','WindowStart_h','WindowEnd_h', ...
           'N_PhaseObs','N_Mice','N_Candidates','N_Transitions','R','MeanPhase_rad','CircularSD_rad', ...
           'MeanRidgePeriod_h','MeanRidgePower_log10','PassN'};

    for g = 1:height(groups)
        gm = TL.Photoperiod_h == groups.Photoperiod_h(g) & ...
             string(TL.BandName) == string(groups.BandName(g)) & ...
             string(TL.TransitionType) == string(groups.TransitionType(g));
        for w = 1:2
            if w == 1
                wStart = -summaryWindowH; wEnd = 0;
                idx = gm & TL.RelativeTime_h >= wStart & TL.RelativeTime_h < wEnd;
            else
                wStart = 0; wEnd = summaryWindowH;
                idx = gm & TL.RelativeTime_h >= wStart & TL.RelativeTime_h <= wEnd;
            end
            theta = TL.RidgePhase_rad(idx);
            theta = theta(isfinite(theta));
            n = numel(theta);
            [R, mu, csd] = circ_summary(theta);
            passN = n >= minN;
            rows(end+1,:) = {groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), winNames(w), wStart, wEnd, ...
                n, count_unique(TL.SignalID(idx)), count_unique(TL.CandidateID(idx)), count_unique(string(TL.CandidateID(idx)) + "|" + string(TL.TransitionType(idx)) + "|" + string(TL.TransitionDay(idx))), ...
                R, mu, csd, mean(TL.RidgePeriod_h(idx),'omitnan'), mean(TL.RidgePower_log10(idx),'omitnan'), passN}; %#ok<AGROW>
        end
    end

    S = rows_to_table(rows, hdr);
end

function D = make_deltaR_summary(PP)
    groups = unique(PP(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    rows = {};
    hdr = {'Photoperiod_h','BandName','TransitionType','R_Pre','R_Post','DeltaR_PostMinusPre', ...
           'MeanPhase_Pre_rad','MeanPhase_Post_rad','N_Pre','N_Post','PassN_Pre','PassN_Post'};
    for g = 1:height(groups)
        idx = PP.Photoperiod_h == groups.Photoperiod_h(g) & ...
              string(PP.BandName) == string(groups.BandName(g)) & ...
              string(PP.TransitionType) == string(groups.TransitionType(g));
        pre = PP(idx & string(PP.Window)=="Pre", :);
        post = PP(idx & string(PP.Window)=="Post", :);
        [rpre,mupre,npre,passpre] = row_vals(pre);
        [rpost,mupost,npost,passpost] = row_vals(post);
        rows(end+1,:) = {groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), ...
            rpre, rpost, rpost-rpre, mupre, mupost, npre, npost, passpre, passpost}; %#ok<AGROW>
    end
    D = rows_to_table(rows, hdr);
end

function S = make_candidate_prepost_summary(TL, summaryWindowH, minN)
    groups = unique(TL(:, {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    winNames = ["Pre"; "Post"];
    rows = {};
    hdr = {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType','Window', ...
           'N_PhaseObs','N_Transitions','R','MeanPhase_rad','CircularSD_rad','MeanRidgePeriod_h','MeanRidgePower_log10','PassN'};
    for g = 1:height(groups)
        gm = string(TL.CandidateID) == string(groups.CandidateID(g)) & ...
             string(TL.TransitionType) == string(groups.TransitionType(g));
        for w = 1:2
            if w == 1
                idx = gm & TL.RelativeTime_h >= -summaryWindowH & TL.RelativeTime_h < 0;
            else
                idx = gm & TL.RelativeTime_h >= 0 & TL.RelativeTime_h <= summaryWindowH;
            end
            theta = TL.RidgePhase_rad(idx); theta = theta(isfinite(theta));
            [R, mu, csd] = circ_summary(theta);
            rows(end+1,:) = {string(groups.CandidateID(g)), string(groups.File(g)), string(groups.SignalID(g)), string(groups.ConditionParsed(g)), ...
                groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), winNames(w), ...
                numel(theta), count_unique(TL.TransitionDay(idx)), R, mu, csd, mean(TL.RidgePeriod_h(idx),'omitnan'), ...
                mean(TL.RidgePower_log10(idx),'omitnan'), numel(theta) >= minN}; %#ok<AGROW>
        end
    end
    S = rows_to_table(rows, hdr);
end

function D = make_candidate_deltaR_summary(CP)
    groups = unique(CP(:, {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    rows = {};
    hdr = {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType','R_Pre','R_Post','DeltaR_PostMinusPre','N_Pre','N_Post','PassN_Pre','PassN_Post'};
    for g = 1:height(groups)
        idx = string(CP.CandidateID) == string(groups.CandidateID(g)) & string(CP.TransitionType) == string(groups.TransitionType(g));
        pre = CP(idx & string(CP.Window)=="Pre", :);
        post = CP(idx & string(CP.Window)=="Post", :);
        [rpre,~,npre,passpre] = row_vals(pre);
        [rpost,~,npost,passpost] = row_vals(post);
        rows(end+1,:) = {string(groups.CandidateID(g)), string(groups.File(g)), string(groups.SignalID(g)), string(groups.ConditionParsed(g)), ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), ...
            rpre, rpost, rpost-rpre, npre, npost, passpre, passpost}; %#ok<AGROW>
    end
    D = rows_to_table(rows, hdr);
end


function Stats = make_resync_primary_stats(CD, minCandidates, nPerm)
    % Candidate-level one-sample directional test: DeltaR > 0 for real transitions.
    hdr = {'FDR_Family','TestName','TestDescription','Photoperiod_h','BandName','TransitionType', ...
           'N_Candidates','N_Mice','MeanDeltaR','MedianDeltaR','SD_DeltaR','Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(CD) || height(CD)==0
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    C = CD(ismember(string(CD.TransitionType), ["DL","LD"]) & isfinite(CD.DeltaR_PostMinusPre) & ...
           logical(CD.PassN_Pre) & logical(CD.PassN_Post), :);
    if isempty(C)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    groups = unique(C(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    for g = 1:height(groups)
        idx = C.Photoperiod_h == groups.Photoperiod_h(g) & ...
              string(C.BandName) == string(groups.BandName(g)) & ...
              string(C.TransitionType) == string(groups.TransitionType(g));
        x = C.DeltaR_PostMinusPre(idx);
        x = x(isfinite(x));
        n = numel(x);
        passN = n >= minCandidates;
        p = NaN;
        if passN
            p = paired_or_onesample_perm_p(x, 'right', nPerm);
        end
        rows(end+1,:) = {"Resync_Primary_DeltaR", ...
            sprintf('DeltaR_gt_0_L%s_%s_%s', safe_num(groups.Photoperiod_h(g)), string(groups.BandName(g)), string(groups.TransitionType(g))), ...
            "One-sample candidate-level test that post-transition phase coherence exceeds pre-transition coherence", ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), ...
            n, count_unique(C.SignalID(idx)), mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'), "right", p, nPerm, passN}; %#ok<AGROW>
    end
    Stats = rows_to_table(rows, hdr);
end

function Stats = make_real_vs_pseudo_stats(CD, minCandidates, nPerm)
    % Candidate-level paired comparison: mean real-transition DeltaR minus
    % mean pseudo-transition DeltaR. Directional test: real > pseudo.
    hdr = {'FDR_Family','TestName','TestDescription','Photoperiod_h','BandName', ...
           'N_Candidates','N_Mice','MeanRealDeltaR','MeanPseudoDeltaR','MeanDifference_RealMinusPseudo', ...
           'MedianDifference_RealMinusPseudo','SD_Difference','Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(CD) || height(CD)==0
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    C = CD(isfinite(CD.DeltaR_PostMinusPre) & logical(CD.PassN_Pre) & logical(CD.PassN_Post), :);
    if isempty(C)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    groups = unique(C(:, {'Photoperiod_h','BandName'}), 'rows', 'stable');
    for g = 1:height(groups)
        G = C(C.Photoperiod_h == groups.Photoperiod_h(g) & string(C.BandName)==string(groups.BandName(g)), :);
        cand = unique(string(G.CandidateID), 'stable');
        realVals = [];
        pseudoVals = [];
        mice = strings(0,1);
        for c = 1:numel(cand)
            Gi = G(string(G.CandidateID)==cand(c), :);
            real = Gi.DeltaR_PostMinusPre(ismember(string(Gi.TransitionType), ["DL","LD"]));
            pseudo = Gi.DeltaR_PostMinusPre(ismember(string(Gi.TransitionType), ["MidLight","MidDark"]));
            if any(isfinite(real)) && any(isfinite(pseudo))
                realVals(end+1,1) = mean(real,'omitnan'); %#ok<AGROW>
                pseudoVals(end+1,1) = mean(pseudo,'omitnan'); %#ok<AGROW>
                mice(end+1,1) = string(Gi.SignalID(1)); %#ok<AGROW>
            end
        end
        diffVals = realVals - pseudoVals;
        n = numel(diffVals);
        passN = n >= minCandidates;
        p = NaN;
        if passN
            p = paired_or_onesample_perm_p(diffVals, 'right', nPerm);
        end
        rows(end+1,:) = {"Resync_RealVsPseudo", ...
            sprintf('Real_gt_Pseudo_L%s_%s', safe_num(groups.Photoperiod_h(g)), string(groups.BandName(g))), ...
            "Paired candidate-level test comparing real LD/DL transition DeltaR against mid-light/mid-dark pseudo-transition DeltaR", ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), ...
            n, count_unique(mice), mean(realVals,'omitnan'), mean(pseudoVals,'omitnan'), mean(diffVals,'omitnan'), ...
            median(diffVals,'omitnan'), std(diffVals,'omitnan'), "right", p, nPerm, passN}; %#ok<AGROW>
    end
    Stats = rows_to_table(rows, hdr);
end

function Stats = make_candidate_metric_prepost_stats(CP, metricName, familyName, tail, minCandidates, nPerm)
    % Candidate-level paired pre/post test for ridge-period or ridge-power metrics.
    hdr = {'FDR_Family','TestName','TestDescription','Metric','Photoperiod_h','BandName','TransitionType', ...
           'N_Candidates','N_Mice','Mean_Pre','Mean_Post','MeanDifference_PostMinusPre', ...
           'MedianDifference_PostMinusPre','SD_Difference','Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(CP) || height(CP)==0 || ~ismember(metricName, CP.Properties.VariableNames)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    C = CP(logical(CP.PassN) & isfinite(CP.(metricName)), :);
    if isempty(C)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    groups = unique(C(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    for g = 1:height(groups)
        G = C(C.Photoperiod_h == groups.Photoperiod_h(g) & ...
              string(C.BandName)==string(groups.BandName(g)) & ...
              string(C.TransitionType)==string(groups.TransitionType(g)), :);
        cand = unique(string(G.CandidateID), 'stable');
        preVals = [];
        postVals = [];
        mice = strings(0,1);
        for c = 1:numel(cand)
            Gi = G(string(G.CandidateID)==cand(c), :);
            pre = Gi.(metricName)(string(Gi.Window)=="Pre");
            post = Gi.(metricName)(string(Gi.Window)=="Post");
            if any(isfinite(pre)) && any(isfinite(post))
                preVals(end+1,1) = mean(pre,'omitnan'); %#ok<AGROW>
                postVals(end+1,1) = mean(post,'omitnan'); %#ok<AGROW>
                mice(end+1,1) = string(Gi.SignalID(1)); %#ok<AGROW>
            end
        end
        diffVals = postVals - preVals;
        n = numel(diffVals);
        passN = n >= minCandidates;
        p = NaN;
        if passN
            p = paired_or_onesample_perm_p(diffVals, tail, nPerm);
        end
        rows(end+1,:) = {string(familyName), ...
            sprintf('%s_PostMinusPre_L%s_%s_%s', metricName, safe_num(groups.Photoperiod_h(g)), string(groups.BandName(g)), string(groups.TransitionType(g))), ...
            "Paired candidate-level pre/post transition test", string(metricName), ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), ...
            n, count_unique(mice), mean(preVals,'omitnan'), mean(postVals,'omitnan'), mean(diffVals,'omitnan'), ...
            median(diffVals,'omitnan'), std(diffVals,'omitnan'), string(tail), p, nPerm, passN}; %#ok<AGROW>
    end
    Stats = rows_to_table(rows, hdr);
end


function Stats = make_projected_primary_stats(CD, minCandidates, nPerm)
    hdr = {'FDR_Family','TestName','TestDescription','Photoperiod_h','BandName','TransitionType', ...
           'N_Candidates','N_Mice','MeanDeltaR','MedianDeltaR','SD_DeltaR','Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(CD) || height(CD)==0
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    C = CD(ismember(string(CD.TransitionType), ["ProjectedDL_LL","ProjectedLD_LL"]) & isfinite(CD.DeltaR_PostMinusPre) & ...
           logical(CD.PassN_Pre) & logical(CD.PassN_Post), :);
    if isempty(C)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    groups = unique(C(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    for g = 1:height(groups)
        idx = C.Photoperiod_h == groups.Photoperiod_h(g) & ...
              string(C.BandName) == string(groups.BandName(g)) & ...
              string(C.TransitionType) == string(groups.TransitionType(g));
        x = C.DeltaR_PostMinusPre(idx);
        x = x(isfinite(x));
        n = numel(x);
        passN = n >= minCandidates;
        p = NaN;
        if passN
            p = paired_or_onesample_perm_p(x, 'right', nPerm);
        end
        rows(end+1,:) = {"LL_Projected_Primary_DeltaR", ...
            sprintf('ProjectedDeltaR_gt_0_L%s_%s_%s', safe_num(groups.Photoperiod_h(g)), string(groups.BandName(g)), string(groups.TransitionType(g))), ...
            "Secondary LL aftereffect test: projected post-event phase coherence exceeds projected pre-event coherence", ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), ...
            n, count_unique(C.SignalID(idx)), mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'), "right", p, nPerm, passN}; %#ok<AGROW>
    end
    Stats = rows_to_table(rows, hdr);
end

function Stats = make_projected_real_vs_pseudo_stats(CD, minCandidates, nPerm)
    hdr = {'FDR_Family','TestName','TestDescription','Photoperiod_h','BandName', ...
           'N_Candidates','N_Mice','MeanProjectedTransitionDeltaR','MeanProjectedPseudoDeltaR','MeanDifference_ProjectedTransitionMinusPseudo', ...
           'MedianDifference_ProjectedTransitionMinusPseudo','SD_Difference','Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(CD) || height(CD)==0
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    C = CD(isfinite(CD.DeltaR_PostMinusPre) & logical(CD.PassN_Pre) & logical(CD.PassN_Post), :);
    if isempty(C)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    groups = unique(C(:, {'Photoperiod_h','BandName'}), 'rows', 'stable');
    for g = 1:height(groups)
        G = C(C.Photoperiod_h == groups.Photoperiod_h(g) & string(C.BandName)==string(groups.BandName(g)), :);
        cand = unique(string(G.CandidateID), 'stable');
        realVals = [];
        pseudoVals = [];
        mice = strings(0,1);
        for c = 1:numel(cand)
            Gi = G(string(G.CandidateID)==cand(c), :);
            real = Gi.DeltaR_PostMinusPre(ismember(string(Gi.TransitionType), ["ProjectedDL_LL","ProjectedLD_LL"]));
            pseudo = Gi.DeltaR_PostMinusPre(ismember(string(Gi.TransitionType), ["ProjectedMidLight_LL","ProjectedMidDark_LL"]));
            if any(isfinite(real)) && any(isfinite(pseudo))
                realVals(end+1,1) = mean(real,'omitnan'); %#ok<AGROW>
                pseudoVals(end+1,1) = mean(pseudo,'omitnan'); %#ok<AGROW>
                mice(end+1,1) = string(Gi.SignalID(1)); %#ok<AGROW>
            end
        end
        diffVals = realVals - pseudoVals;
        n = numel(diffVals);
        passN = n >= minCandidates;
        p = NaN;
        if passN
            p = paired_or_onesample_perm_p(diffVals, 'right', nPerm);
        end
        rows(end+1,:) = {"LL_Projected_TransitionVsPseudo", ...
            sprintf('ProjectedTransition_gt_Pseudo_L%s_%s', safe_num(groups.Photoperiod_h(g)), string(groups.BandName(g))), ...
            "Secondary LL aftereffect test comparing projected former-transition DeltaR against projected pseudo-transition DeltaR", ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), ...
            n, count_unique(mice), mean(realVals,'omitnan'), mean(pseudoVals,'omitnan'), mean(diffVals,'omitnan'), ...
            median(diffVals,'omitnan'), std(diffVals,'omitnan'), "right", p, nPerm, passN}; %#ok<AGROW>
    end
    Stats = rows_to_table(rows, hdr);
end

function S = make_candidate_day_prepost_summary(TL, summaryWindowH, minN)
    groups = unique(TL(:, {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType','TransitionDay'}), 'rows', 'stable');
    winNames = ["Pre"; "Post"];
    rows = {};
    hdr = {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType','TransitionDay','Window', ...
           'N_PhaseObs','R','MeanPhase_rad','CircularSD_rad','MeanRidgePeriod_h','MeanRidgePower_log10','PassN'};
    for g = 1:height(groups)
        gm = string(TL.CandidateID) == string(groups.CandidateID(g)) & ...
             string(TL.TransitionType) == string(groups.TransitionType(g)) & ...
             TL.TransitionDay == groups.TransitionDay(g);
        for w = 1:2
            if w == 1
                idx = gm & TL.RelativeTime_h >= -summaryWindowH & TL.RelativeTime_h < 0;
            else
                idx = gm & TL.RelativeTime_h >= 0 & TL.RelativeTime_h <= summaryWindowH;
            end
            theta = TL.RidgePhase_rad(idx); theta = theta(isfinite(theta));
            [R, mu, csd] = circ_summary(theta);
            rows(end+1,:) = {string(groups.CandidateID(g)), string(groups.File(g)), string(groups.SignalID(g)), string(groups.ConditionParsed(g)), ...
                groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), groups.TransitionDay(g), winNames(w), ...
                numel(theta), R, mu, csd, mean(TL.RidgePeriod_h(idx),'omitnan'), mean(TL.RidgePower_log10(idx),'omitnan'), numel(theta) >= minN}; %#ok<AGROW>
        end
    end
    S = rows_to_table(rows, hdr);
end

function D = make_candidate_day_deltaR_summary(CP)
    groups = unique(CP(:, {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType','TransitionDay'}), 'rows', 'stable');
    rows = {};
    hdr = {'CandidateID','File','SignalID','ConditionParsed','Photoperiod_h','BandName','TransitionType','TransitionDay','R_Pre','R_Post','DeltaR_PostMinusPre','N_Pre','N_Post','PassN_Pre','PassN_Post'};
    for g = 1:height(groups)
        idx = string(CP.CandidateID) == string(groups.CandidateID(g)) & ...
              string(CP.TransitionType) == string(groups.TransitionType(g)) & ...
              CP.TransitionDay == groups.TransitionDay(g);
        pre = CP(idx & string(CP.Window)=="Pre", :);
        post = CP(idx & string(CP.Window)=="Post", :);
        [rpre,~,npre,passpre] = row_vals(pre);
        [rpost,~,npost,passpost] = row_vals(post);
        rows(end+1,:) = {string(groups.CandidateID(g)), string(groups.File(g)), string(groups.SignalID(g)), string(groups.ConditionParsed(g)), ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), groups.TransitionDay(g), ...
            rpre, rpost, rpost-rpre, npre, npost, passpre, passpost}; %#ok<AGROW>
    end
    D = rows_to_table(rows, hdr);
end

function Stats = make_ll_projected_day_decay_stats(DayDelta, minCandidates, nPerm)
    hdr = {'FDR_Family','TestName','TestDescription','Photoperiod_h','BandName','TransitionType', ...
           'N_DayCandidateObservations','N_Candidates','N_Mice','SpearmanRho_DayVsDeltaR','Slope_DeltaRPerDay', ...
           'Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(DayDelta) || height(DayDelta)==0
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    C = DayDelta(isfinite(DayDelta.DeltaR_PostMinusPre) & logical(DayDelta.PassN_Pre) & logical(DayDelta.PassN_Post), :);
    C = C(ismember(string(C.TransitionType), ["ProjectedDL_LL","ProjectedLD_LL"]), :);
    if isempty(C)
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    groups = unique(C(:, {'Photoperiod_h','BandName','TransitionType'}), 'rows', 'stable');
    for g = 1:height(groups)
        idx = C.Photoperiod_h == groups.Photoperiod_h(g) & ...
              string(C.BandName)==string(groups.BandName(g)) & ...
              string(C.TransitionType)==string(groups.TransitionType(g));
        x = double(C.TransitionDay(idx));
        y = double(C.DeltaR_PostMinusPre(idx));
        valid = isfinite(x) & isfinite(y);
        x = x(valid); y = y(valid);
        n = numel(y);
        passN = n >= minCandidates && numel(unique(x)) >= 2;
        rho = NaN; slope = NaN; p = NaN;
        if n >= 2
            rho = corr(x(:), y(:), 'Type','Spearman', 'Rows','complete');
            pp = polyfit(x(:), y(:), 1); slope = pp(1);
        end
        if passN
            p = spearman_perm_p(x, y, 'two-sided', nPerm);
        end
        rows(end+1,:) = {"LL_Projected_DayDecay", ...
            sprintf('LLProjected_DayDecay_L%s_%s_%s', safe_num(groups.Photoperiod_h(g)), string(groups.BandName(g)), string(groups.TransitionType(g))), ...
            "Secondary LL aftereffect test: projected DeltaR association with day in LL", ...
            groups.Photoperiod_h(g), string(groups.BandName(g)), string(groups.TransitionType(g)), ...
            n, count_unique(C.CandidateID(idx)), count_unique(C.SignalID(idx)), rho, slope, "two-sided", p, nPerm, passN}; %#ok<AGROW>
    end
    Stats = rows_to_table(rows, hdr);
end

function Stats = make_l22_vs_ll_projected_stats(CDReal, CDLL, refPhoto, minCandidates, nPerm)
    hdr = {'FDR_Family','TestName','TestDescription','BandName','TransitionPair','ReferencePhotoperiod_h', ...
           'N_L22','N_LL','MeanDeltaR_L22','MeanDeltaR_LLProjected','MeanDifference_L22MinusLLProjected', ...
           'MedianDifference_NA','Tail','PValue_raw','PermutationN','PassN'};
    rows = {};
    if isempty(CDReal) || isempty(CDLL) || height(CDReal)==0 || height(CDLL)==0
        Stats = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr); return;
    end
    pairs = {"DL", "ProjectedDL_LL"; "LD", "ProjectedLD_LL"};
    bands = unique([string(CDReal.BandName); string(CDLL.BandName)], 'stable');
    for b = 1:numel(bands)
        for p = 1:size(pairs,1)
            realType = pairs{p,1}; llType = pairs{p,2};
            x = CDReal.DeltaR_PostMinusPre(CDReal.Photoperiod_h == refPhoto & string(CDReal.BandName)==bands(b) & string(CDReal.TransitionType)==realType & isfinite(CDReal.DeltaR_PostMinusPre) & logical(CDReal.PassN_Pre) & logical(CDReal.PassN_Post));
            y = CDLL.DeltaR_PostMinusPre(string(CDLL.BandName)==bands(b) & string(CDLL.TransitionType)==llType & isfinite(CDLL.DeltaR_PostMinusPre) & logical(CDLL.PassN_Pre) & logical(CDLL.PassN_Post));
            nX = numel(x); nY = numel(y);
            passN = nX >= minCandidates && nY >= minCandidates;
            pval = NaN;
            if passN
                pval = twosample_perm_p(x, y, 'two-sided', nPerm);
            end
            rows(end+1,:) = {"L22_vs_LLProjected_Aftereffect", ...
                sprintf('L22Real_vs_LLProjected_%s_%s', string(realType), bands(b)), ...
                "Secondary aftereffect comparison: final real L22 transition DeltaR versus projected former-transition DeltaR in LL", ...
                bands(b), sprintf('%s_vs_%s', string(realType), string(llType)), refPhoto, ...
                nX, nY, mean(x,'omitnan'), mean(y,'omitnan'), mean(x,'omitnan')-mean(y,'omitnan'), ...
                NaN, "two-sided", pval, nPerm, passN}; %#ok<AGROW>
        end
    end
    Stats = rows_to_table(rows, hdr);
end

function p = twosample_perm_p(x, y, tail, nPerm)
    x = double(x(:)); y = double(y(:));
    x = x(isfinite(x)); y = y(isfinite(y));
    if isempty(x) || isempty(y), p = NaN; return; end
    obs = mean(x,'omitnan') - mean(y,'omitnan');
    pool = [x; y]; nx = numel(x); n = numel(pool);
    if nPerm <= 0 || ~isfinite(nPerm), nPerm = 10000; end
    stats = nan(nPerm,1);
    for i = 1:nPerm
        ord = randperm(n);
        xi = pool(ord(1:nx)); yi = pool(ord(nx+1:end));
        stats(i) = mean(xi,'omitnan') - mean(yi,'omitnan');
    end
    switch lower(string(tail))
        case "right"
            p = (1 + sum(stats >= obs)) / (nPerm + 1);
        case "left"
            p = (1 + sum(stats <= obs)) / (nPerm + 1);
        otherwise
            p = (1 + sum(abs(stats) >= abs(obs))) / (nPerm + 1);
    end
end

function p = spearman_perm_p(x, y, tail, nPerm)
    x = double(x(:)); y = double(y(:));
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);
    if numel(y) < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
        p = NaN; return;
    end
    obs = corr(x, y, 'Type','Spearman', 'Rows','complete');
    if nPerm <= 0 || ~isfinite(nPerm), nPerm = 10000; end
    stats = nan(nPerm,1);
    for i = 1:nPerm
        yp = y(randperm(numel(y)));
        stats(i) = corr(x, yp, 'Type','Spearman', 'Rows','complete');
    end
    switch lower(string(tail))
        case "right"
            p = (1 + sum(stats >= obs)) / (nPerm + 1);
        case "left"
            p = (1 + sum(stats <= obs)) / (nPerm + 1);
        otherwise
            p = (1 + sum(abs(stats) >= abs(obs))) / (nPerm + 1);
    end
end

function p = paired_or_onesample_perm_p(x, tail, nPerm)
    % Sign-flip permutation test for paired differences or one-sample values.
    % tail: 'right', 'left', or 'two-sided'. Statistic is mean(x).
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        p = NaN; return;
    end
    obs = mean(x, 'omitnan');
    if all(x == 0)
        p = 1; return;
    end
    n = numel(x);
    if nargin < 3 || isempty(nPerm) || ~isfinite(nPerm) || nPerm <= 0
        nPerm = 10000;
    end
    % Exact enumeration is feasible for small n.
    if n <= 15
        combos = dec2bin(0:(2^n-1)) - '0';
        signs = 2*combos - 1;
        permStats = mean(signs .* reshape(x,1,[]), 2);
    else
        signs = 2*(rand(nPerm,n) > 0.5) - 1;
        permStats = mean(signs .* reshape(x,1,[]), 2);
    end
    switch lower(string(tail))
        case "right"
            p = (1 + sum(permStats >= obs)) / (numel(permStats) + 1);
        case "left"
            p = (1 + sum(permStats <= obs)) / (numel(permStats) + 1);
        otherwise
            p = (1 + sum(abs(permStats) >= abs(obs))) / (numel(permStats) + 1);
    end
end

function Settings = make_settings_table(scriptName, scriptVersion, timestamp, mapPath, handoffDir, primaryPhase, requireValidFlag, include12h, urBands, periWindow, binWidth, summaryWindow, minBin, minSummary, alphaFDR, nPermStats, minCandidatesForTest, adaptiveWindows, windowFraction, excludeInitialDays, minDayStable, llPhotoValue)
    Settings = table();
    Settings.ScriptName = string(scriptName);
    Settings.ScriptVersion = string(scriptVersion);
    Settings.Timestamp = string(timestamp);
    Settings.ValidationMap = string(mapPath);
    Settings.HandoffFolder = string(handoffDir);
    Settings.PrimaryPhase = string(primaryPhase);
    Settings.RequireValidFlag = logical(requireValidFlag);
    Settings.IncludeHarmonicSensitive12h = logical(include12h);
    Settings.URBands = join(string(urBands), '; ');
    Settings.PeriTransitionWindow_h = periWindow;
    Settings.BinWidth_h = binWidth;
    Settings.PrePostSummaryWindow_h = summaryWindow;
    Settings.MinPointsPerBin = minBin;
    Settings.MinPointsSummary = minSummary;
    Settings.AlphaFDR = alphaFDR;
    Settings.NPermStats = nPermStats;
    Settings.MinCandidatesForTest = minCandidatesForTest;
    Settings.AdaptiveTransitionWindows = logical(adaptiveWindows);
    Settings.TransitionWindowFraction = windowFraction;
    Settings.ExcludeInitialDaysFromResync = logical(excludeInitialDays);
    Settings.MinDayForStableResync = minDayStable;
    Settings.LLPhotoperiodValue = llPhotoValue;
    Settings.LLHandling = "Photoperiod_h >= LLPhotoperiodValue is excluded from LD/DL transition inference and reported in LL_NoTransitionSummary";
    Settings.FDRMethod = "Benjamini-Hochberg applied within predefined inference families";
    Settings.PrimaryInference = "Candidate-level DeltaR > 0 for real DL/LD transitions";
end

function make_resync_figures(Bin, DeltaR, CandDelta, figDir, dpi, ext)
    ensure_dir(figDir);

    % Figure 1: phase coherence curves, one file per photoperiod-band.
    if ~isempty(Bin) && height(Bin) > 0
        keep = Bin.PassN & isfinite(Bin.R);
        B = Bin(keep,:);
        combos = unique(B(:, {'Photoperiod_h','BandName'}), 'rows', 'stable');
        for i = 1:height(combos)
            idx = B.Photoperiod_h == combos.Photoperiod_h(i) & string(B.BandName) == string(combos.BandName(i));
            Bi = B(idx,:);
            if isempty(Bi), continue; end
            f = figure('Color','w','Position',[100 100 950 520]);
            hold on;
            ttypes = ["DL","LD","MidLight","MidDark"];
            for t = 1:numel(ttypes)
                Bt = Bi(string(Bi.TransitionType)==ttypes(t),:);
                if isempty(Bt), continue; end
                Bt = sortrows(Bt, 'RelBinCenter_h');
                plot(Bt.RelBinCenter_h, Bt.R, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', char(ttypes(t)));
            end
            xline(0, ':', 'Transition', 'LineWidth', 1.2);
            xlabel('Time relative to event (h)', 'FontWeight','bold');
            ylabel('Phase coherence, R', 'FontWeight','bold');
            title(sprintf('Ridge-following phase coherence | L%.3g | %s', combos.Photoperiod_h(i), string(combos.BandName(i))), 'Interpreter','none');
            legend('Location','eastoutside', 'Interpreter','none');
            ylim([0 1]); box off; set(gca,'TickDir','out','FontName','Times New Roman');
            out = fullfile(figDir, sprintf('PhaseCoherenceCurve_L%s_%s%s', safe_num(combos.Photoperiod_h(i)), sanitise_filename(string(combos.BandName(i))), ext));
            exportgraphics(f, out, 'Resolution', dpi);
            close(f);
        end
    end

    % Figure 2: DeltaR grouped by photoperiod-band.
    if ~isempty(DeltaR) && height(DeltaR) > 0
        combos = unique(DeltaR(:, {'Photoperiod_h','BandName'}), 'rows', 'stable');
        for i = 1:height(combos)
            idx = DeltaR.Photoperiod_h == combos.Photoperiod_h(i) & string(DeltaR.BandName) == string(combos.BandName(i));
            D = DeltaR(idx,:);
            if isempty(D), continue; end
            order = ["DL","LD","MidLight","MidDark"];
            vals = nan(numel(order),1);
            for j = 1:numel(order)
                k = find(string(D.TransitionType)==order(j),1);
                if ~isempty(k), vals(j)=D.DeltaR_PostMinusPre(k); end
            end
            f = figure('Color','w','Position',[100 100 720 470]);
            bar(vals);
            set(gca,'XTick',1:numel(order),'XTickLabel',cellstr(order),'TickDir','out','FontName','Times New Roman');
            yline(0, ':', 'LineWidth', 1.2);
            ylabel('\DeltaR (post - pre)', 'FontWeight','bold');
            title(sprintf('Transition phase-resynchronisation | L%.3g | %s', combos.Photoperiod_h(i), string(combos.BandName(i))), 'Interpreter','none');
            box off;
            out = fullfile(figDir, sprintf('DeltaR_L%s_%s%s', safe_num(combos.Photoperiod_h(i)), sanitise_filename(string(combos.BandName(i))), ext));
            exportgraphics(f, out, 'Resolution', dpi);
            close(f);
        end
    end

    % Figure 3: candidate-level DeltaR distributions across real vs pseudo events.
    if ~isempty(CandDelta) && height(CandDelta) > 0
        C = CandDelta(isfinite(CandDelta.DeltaR_PostMinusPre), :);
        if ~isempty(C)
            f = figure('Color','w','Position',[100 100 900 520]);
            cats = categorical(string(C.TransitionType), ["DL","LD","MidLight","MidDark"]);
            boxchart(cats, C.DeltaR_PostMinusPre);
            yline(0, ':', 'LineWidth', 1.2);
            ylabel('Candidate-level \DeltaR (post - pre)', 'FontWeight','bold');
            xlabel('Event type', 'FontWeight','bold');
            title('Candidate-level ridge-phase resynchronisation', 'Interpreter','none');
            box off; set(gca,'TickDir','out','FontName','Times New Roman');
            exportgraphics(f, fullfile(figDir, ['CandidateLevel_DeltaR_ByEvent' ext]), 'Resolution', dpi);
            close(f);
        end
    end
end


function make_ll_projected_figures(Bin, DeltaR, CandDelta, figDir, dpi, ext)
    ensure_dir(figDir);
    outDir = fullfile(figDir, 'LL_Projected_Aftereffect');
    ensure_dir(outDir);

    if ~isempty(Bin) && height(Bin) > 0
        keep = Bin.PassN & isfinite(Bin.R);
        B = Bin(keep,:);
        combos = unique(B(:, {'Photoperiod_h','BandName'}), 'rows', 'stable');
        order = ["ProjectedDL_LL","ProjectedLD_LL","ProjectedMidLight_LL","ProjectedMidDark_LL"];
        for i = 1:height(combos)
            idx = B.Photoperiod_h == combos.Photoperiod_h(i) & string(B.BandName) == string(combos.BandName(i));
            Bi = B(idx,:);
            if isempty(Bi), continue; end
            f = figure('Color','w','Position',[100 100 1000 540]);
            hold on;
            for t = 1:numel(order)
                Bt = Bi(string(Bi.TransitionType)==order(t),:);
                if isempty(Bt), continue; end
                Bt = sortrows(Bt, 'RelBinCenter_h');
                plot(Bt.RelBinCenter_h, Bt.R, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', char(order(t)));
            end
            xline(0, ':', 'Projected event', 'LineWidth', 1.2);
            xlabel('Time relative to projected former event (h)', 'FontWeight','bold');
            ylabel('Phase coherence, R', 'FontWeight','bold');
            title(sprintf('LL projected ridge-phase coherence | %s', string(combos.BandName(i))), 'Interpreter','none');
            legend('Location','eastoutside', 'Interpreter','none');
            ylim([0 1]); box off; set(gca,'TickDir','out','FontName','Times New Roman');
            out = fullfile(outDir, sprintf('LLProjected_PhaseCoherence_%s%s', sanitise_filename(string(combos.BandName(i))), ext));
            exportgraphics(f, out, 'Resolution', dpi);
            close(f);
        end
    end

    if ~isempty(DeltaR) && height(DeltaR) > 0
        combos = unique(DeltaR(:, {'Photoperiod_h','BandName'}), 'rows', 'stable');
        order = ["ProjectedDL_LL","ProjectedLD_LL","ProjectedMidLight_LL","ProjectedMidDark_LL"];
        for i = 1:height(combos)
            idx = DeltaR.Photoperiod_h == combos.Photoperiod_h(i) & string(DeltaR.BandName) == string(combos.BandName(i));
            D = DeltaR(idx,:);
            vals = nan(numel(order),1);
            for j = 1:numel(order)
                k = find(string(D.TransitionType)==order(j),1);
                if ~isempty(k), vals(j)=D.DeltaR_PostMinusPre(k); end
            end
            f = figure('Color','w','Position',[100 100 820 470]);
            bar(vals);
            set(gca,'XTick',1:numel(order),'XTickLabel',cellstr(order),'TickDir','out','FontName','Times New Roman');
            xtickangle(25);
            yline(0, ':', 'LineWidth', 1.2);
            ylabel('\DeltaR (post - pre)', 'FontWeight','bold');
            title(sprintf('LL projected aftereffect | %s', string(combos.BandName(i))), 'Interpreter','none');
            box off;
            out = fullfile(outDir, sprintf('LLProjected_DeltaR_%s%s', sanitise_filename(string(combos.BandName(i))), ext));
            exportgraphics(f, out, 'Resolution', dpi);
            close(f);
        end
    end

    if ~isempty(CandDelta) && height(CandDelta) > 0
        C = CandDelta(isfinite(CandDelta.DeltaR_PostMinusPre), :);
        if ~isempty(C)
            f = figure('Color','w','Position',[100 100 980 540]);
            cats = categorical(string(C.TransitionType), ["ProjectedDL_LL","ProjectedLD_LL","ProjectedMidLight_LL","ProjectedMidDark_LL"]);
            boxchart(cats, C.DeltaR_PostMinusPre);
            yline(0, ':', 'LineWidth', 1.2);
            ylabel('Candidate-level \DeltaR (post - pre)', 'FontWeight','bold');
            xlabel('Projected former event type', 'FontWeight','bold');
            title('LL candidate-level projected ridge-phase aftereffect', 'Interpreter','none');
            box off; set(gca,'TickDir','out','FontName','Times New Roman');
            xtickangle(25);
            exportgraphics(f, fullfile(outDir, ['LLProjected_CandidateLevel_DeltaR_ByEvent' ext]), 'Resolution', dpi);
            close(f);
        end
    end
end

%% ----------------------------- SMALL HELPERS ----------------------------
function [R, mu, csd] = circ_summary(theta)
    theta = theta(:);
    theta = theta(isfinite(theta));
    if isempty(theta)
        R = NaN; mu = NaN; csd = NaN; return;
    end
    z = mean(exp(1i*theta));
    R = abs(z);
    mu = angle(z);
    if R > 0
        csd = sqrt(max(0, -2*log(R)));
    else
        csd = Inf;
    end
end

function [R, mu, n, passN] = row_vals(T)
    if isempty(T)
        R = NaN; mu = NaN; n = 0; passN = false; return;
    end
    R = T.R(1);
    mu = T.MeanPhase_rad(1);
    n = T.N_PhaseObs(1);
    passN = logical(T.PassN(1));
end


function T = rows_to_table(rows, hdr)
    if isempty(rows)
        T = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr);
    else
        T = cell2table(rows, 'VariableNames', hdr);
    end
end

function n = count_unique(x)
    if isempty(x), n = 0; else, n = numel(unique(string(x))); end
end

function tf = to_logical(x)
    if islogical(x)
        tf = x(:);
    elseif isnumeric(x)
        tf = x(:) ~= 0;
    else
        s = lower(strtrim(string(x(:))));
        tf = ismember(s, ["true","1","yes","y","pass","passed"]);
    end
end

function v = robust_scalar_mode(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        v = NaN; return;
    end
    try
        v = mode(x);
    catch
        v = median(x, 'omitnan');
    end
end

function val = get_optional_num(T, name, r)
    if ismember(name, T.Properties.VariableNames)
        val = double(T.(name)(r));
    else
        val = NaN;
    end
end

function val = get_optional_str(T, name, r)
    if ismember(name, T.Properties.VariableNames)
        val = string(T.(name)(r));
    else
        val = string(missing);
    end
end

function val = get_optional_logical(T, name, r)
    if ismember(name, T.Properties.VariableNames)
        val = logical(T.(name)(r));
    else
        val = false;
    end
end

function T = vertcat_compatible(A, B)
    if isempty(A) || width(A)==0
        T = B; return;
    end
    if isempty(B) || width(B)==0
        T = A; return;
    end
    allVars = unique([A.Properties.VariableNames, B.Properties.VariableNames], 'stable');
    A = add_missing_vars(A, allVars);
    B = add_missing_vars(B, allVars);
    T = [A(:,allVars); B(:,allVars)];
end

function T = add_missing_vars(T, allVars)
    for i = 1:numel(allVars)
        v = allVars{i};
        if ~ismember(v, T.Properties.VariableNames)
            T.(v) = repmat(missing, height(T), 1);
        end
    end
end

function safe_writetable_xlsx(T, xlsxPath, sheetName, outRoot, maxRows)
    if isempty(T)
        writecell({'No rows.'}, xlsxPath, 'Sheet', sheetName);
        return;
    end
    if height(T) <= maxRows
        writetable(T, xlsxPath, 'Sheet', sheetName);
    else
        csvPath = fullfile(outRoot, [sheetName '.csv']);
        writetable(T, csvPath);
        writecell({sprintf('Table exceeded Excel sheet limit. Full table written to CSV: %s', csvPath); ...
                   sprintf('Rows: %d', height(T))}, xlsxPath, 'Sheet', sheetName);
    end
end

function ensure_dir(p)
    if ~exist(p,'dir'), mkdir(p); end
end

function delete_if_exists(p)
    if isfile(p), delete(p); end
end

function log_line(fid, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    fprintf('%s\n', msg);
    if fid > 0
        fprintf(fid, '[%s] %s\n', datestr(now,31), msg);
    end
end

function fclose_if_open(fid)
    if ~isempty(fid) && fid > 0
        fclose(fid);
    end
end

function s = sanitise_filename(s)
    s = char(string(s));
    s = regexprep(s, '[<>:"/\\|?*]', '_');
    s = regexprep(s, '\s+', '_');
    s = regexprep(s, '_+', '_');
    if isempty(s), s = 'unnamed'; end
end

function s = safe_num(x)
    s = sprintf('%.3g', x);
    s = strrep(s, '.', 'p');
    s = strrep(s, '-', 'm');
end

