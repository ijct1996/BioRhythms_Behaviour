function cfg = core_defaults()
%CORE_DEFAULTS Locked Core v1 analysis defaults.
%
%   cfg = core_defaults()

    cfg.version = '1.0';
    cfg.matlabTarget = 'R2025b';

    %% Sampling
    cfg.samplingHours = 10 / 60;          % 10 min → 0.167 h
    cfg.samplingMinutes = 10;

    %% Photoperiod / cohort scope (Core v1)
    cfg.coreCohorts = {'C57_LP', 'NR2B_LP', 'NR2B_LD_DD'};
    cfg.photoperiodStepHours = 2;
    cfg.llHours = 24;                     % L24 = constant light
    cfg.ldHours = 0;                      % L0 = constant darkness
    cfg.includeProjectedDark = false;     % Extended UR only

    %% Harmonic subtraction
    cfg.hsub.minPeriodsMinutes = [360, 60];
    cfg.hsub.defaultResidual = 'SEL_P360'; % Selective Min360 when AnchorOK
    cfg.hsub.anchorBandHours = [22, 28];
    cfg.hsub.anchorStepHours = 0.01;
    cfg.hsub.blockLenHours = 8;
    cfg.hsub.nSurrogates = 300;
    cfg.hsub.alphaAnchor = 0.05;
    cfg.hsub.missingFracThreshold = 0.05;
    cfg.hsub.maxGapMinutes = 30;
    cfg.hsub.minCyclesForAnchor = 3.5;
    cfg.hsub.minDeltaR2 = 0.005;
    cfg.hsub.useEdgeGuard = true;
    cfg.hsub.edgeMarginHours = 0.10;
    cfg.hsub.subtractFundamentalInSelective = true;

    %% Wavelet
    cfg.wavelet.waveletName = 'amor';
    cfg.wavelet.periodLimitsMinutes = [60, 1590];
    cfg.wavelet.scalogramYTicksHours = 0:4:26;
    cfg.wavelet.topNPeaks = 5;
    cfg.wavelet.colormap = 'jet';         % NON-NEGOTIABLE for scalograms

    %% Period / co-expression bands (hours)
    cfg.bands.CR_20_28 = [20, 28];
    cfg.bands.UR_1_3   = [1, 3];
    cfg.bands.UR_3_6   = [3, 6];
    cfg.bands.UR_6_12  = [6, 12];
    cfg.bands.UR_12_18 = [12, 18];
    cfg.bands.coexpression = {'CR_20_28', 'UR_1_3', 'UR_3_6', 'UR_6_12', 'UR_12_18'};

    %% Handoff
    cfg.handoff.version = '1.0';
    cfg.handoff.indexName = 'CoreHandoff_Index.xlsx';
    cfg.handoff.summaryPrefix = 'CoreSummary__';

    %% Plot export
    cfg.plot.development.dpi = 96;
    cfg.plot.development.format = 'png';
    cfg.plot.publication.dpi = 600;
    cfg.plot.publication.format = 'jpeg';
    cfg.plot.scalogram.development = struct('dpi', 96, 'format', 'png', 'colormap', 'jet');
    cfg.plot.scalogram.publication = struct('dpi', 600, 'format', 'jpeg', 'colormap', 'jet');
    cfg.plot.hsubIndividual = struct('dpi', 150, 'format', 'png', 'colormap', 'jet');

    %% HSub figure export (Script 1)
    cfg.hsub.scalogramArm = 'Min360';           % SEL_P360 selective removed/residual
    cfg.hsub.scalogramLabel = 'SEL_P360';
    cfg.hsub.saveScalograms = true;
    cfg.hsub.saveIndividualScalogramsDefault = false;

    %% Script 3 — HSub-validated ultradian (residual CWT; AnchorOK only)
    cfg.hsubValidation.harmonicToleranceFrac = 0.05;   % ±5% of P0/k
    cfg.hsubValidation.maxHarmonicOrder = 6;
    cfg.hsubValidation.urPeriodRangeHours = [1, 18];
    cfg.hsubValidation.topNPeaks = 5;
    cfg.hsubValidation.crSource = 'raw';               % CR_20_28 from RAW; UR from residual

    %% Excel column conventions
    cfg.io.timeColumnCandidates = {'Time (hr)', 'Time', 'Time_hr'};
    cfg.io.lightColumnCandidates = {'Light duration (h)', 'Light duration (h) '};
end
