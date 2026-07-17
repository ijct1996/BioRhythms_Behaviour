function cfg = extended_defaults()
%EXTENDED_DEFAULTS Extended Scripts 4–7 defaults (development-first).
%
%   cfg = extended_defaults()
%
%   Core Scripts 1–3 remain frozen. Extended adds ridge validation,
%   transition resync, LME/FDR, and phase/profile outputs only.

    cfg.version = '1.0-dev';
    cfg.matlabTarget = 'R2025b';
    cfg.phase = 'Scripts4-7';

    %% Plot mode — development until user requests publication pass
    cfg.plotMode = 'development';
    cfg.plot.generateFigures = true;
    cfg.plot.generateLegacyFigures = false;  % Kent E/F/G still 600 JPEG if enabled
    cfg.plot.saveTiff = false;

    %% Sampling (same as Core)
    cfg.samplingHours = 10 / 60;
    cfg.samplingMinutes = 10;

    %% Script numbering (Extended 4–7)
    cfg.scripts = struct();
    cfg.scripts(4).name = 'Ridge handoff + CarryForward validation';
    cfg.scripts(4).entry = 'run_extended_script4_ridge_validation';
    cfg.scripts(5).name = 'Transition resync, FDR, LL projected';
    cfg.scripts(5).entry = 'run_extended_script5_transition_resync';
    cfg.scripts(6).name = 'Across-photoperiod LME / FDR (validated raw UR)';
    cfg.scripts(6).entry = 'run_extended_script6_across_lme';
    cfg.scripts(7).name = 'Phase events + publication profiles (dev tables/figs)';
    cfg.scripts(7).entry = 'run_extended_script7_phase_profiles';

    %% HSub / CarryForward (Script 4 gate)
    cfg.hsub.defaultResidual = 'SEL_P360';
    cfg.hsub.primaryMode = "SEL_P360";
    cfg.hsub.secondaryModes = ["SEL_P60"];
    cfg.hsub.sensitivityModes = ["FL_P360", "FL_P60"];
    cfg.hsub.residualArms = struct('SEL_P360', 'Min360', 'SEL_P60', 'Min60', ...
        'FL_P360', 'Min360', 'FL_P60', 'Min60');
    cfg.hsub.residualTypes = struct('SEL_P360', 'Selective', 'SEL_P60', 'Selective', ...
        'FL_P360', 'FullLadder', 'FL_P60', 'FullLadder');
    cfg.carryForward.periodToleranceFrac = 0.15;
    cfg.carryForward.primaryPhase = "All";
    cfg.carryForward.minRawRidgeCoverage = 0.50;
    cfg.carryForward.minRawCOIValidFrac = 0.50;
    cfg.carryForward.minHSubRidgeCoverage = 0.50;
    cfg.carryForward.minHSubCOIValidFrac = 0.50;
    cfg.carryForward.requireRawPassQC = true;
    cfg.carryForward.requireHSubPassQC = true;

    %% Bands — Kent split; UR_1_3 and UR_3_6 are co-primary for transition story
    cfg.bands.CR = [20, 28];
    cfg.bands.UR = [
        1, 3;
        3, 6;
        6, 9;
        9, 12;
        12, 18];
    cfg.bands.UR_names = ["UR_1_3", "UR_3_6", "UR_6_9", "UR_9_12", "UR_12_18"];
    cfg.bands.primaryUR = ["UR_1_3", "UR_3_6"];
    cfg.bands.harmonicSensitive = ["UR_12_18"];
    cfg.bands.allNames = ["CR_20_28", cfg.bands.UR_names];

    %% Ridge handoff QC (Script 4)
    cfg.ridgeHandoff.minRidgeCoverage = 0.50;
    cfg.ridgeHandoff.minCOIValidFrac = 0.50;

    %% Ridge-phase resync (Script 5)
    cfg.ridge.periWindowH = 6.0;
    cfg.ridge.binWidthH = 0.5;
    cfg.ridge.summaryWindowH = 2.0;
    cfg.ridge.minPointsPerBin = 10;
    cfg.ridge.minPointsSummary = 10;
    cfg.ridge.requireValidFlag = true;
    cfg.ridge.includeHarmonicSensitive12h = true;
    cfg.ridge.adaptiveTransitionWindows = true;
    cfg.ridge.transitionWindowFraction = 0.45;
    cfg.ridge.excludeInitialDays = false;
    cfg.ridge.minDayForStableResync = 0.0;
    cfg.ridge.runStableDaysSensitivity = true;
    cfg.ridge.stableDaysMin = 2.0;
    cfg.ridge.primaryBandsOnly = true;  % gradient figures focus on UR_1_3 + UR_3_6

    %% LL projected dark
    cfg.ll.photoperiodValue = 24.0;
    cfg.ll.projectedReferencePhotoperiodH = 22.0;
    cfg.ll.projectedDarkZT = [22, 24];
    cfg.ll.doProjectedAftereffect = true;
    cfg.ll.runProjectedDayDecay = true;
    cfg.ll.runL22VsProjectedComparison = true;

    %% Inferential statistics
    cfg.stats.alphaFdr = 0.05;
    cfg.stats.nPerm = 10000;
    cfg.stats.minCandidatesForTest = 3;

    %% Plot export — filled by extended_apply_plot_cfg
    cfg.plot.saveDpi = 96;
    cfg.plot.figExt = '.png';

    cfg = extended_apply_plot_cfg(cfg);

    %% Parallel claims (methods lock)
    cfg.parallelClaims.coreScript3 = ...
        'Residual-CWT UR (AnchorOK, ±5% harmonic exclusion) — Core METHODS.md';
    cfg.parallelClaims.extended = ...
        'CarryForward Raw ridge UR (±15% SEL_P360) — Extended METHODS_EXTENDED.md';
end
