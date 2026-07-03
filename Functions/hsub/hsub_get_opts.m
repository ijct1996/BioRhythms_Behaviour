function opts = hsub_get_opts()
%HSUB_GET_OPTS Harmonic subtraction parameters (from core_defaults + QC thresholds).

    cfg = core_defaults();
    opts = cfg.hsub;
    opts.minPeriodsToRun = cfg.hsub.minPeriodsMinutes;
    opts.peakSearchFrac = 0.06;
    opts.ratioTol = 0.08;
    opts.Rcrit = 0.55;
    opts.rhoCrit = 0.45;
    opts.ampFracCritShort = 0.04;
    opts.stepHours = 24;
    opts.minSamplesForAnalysis = 250;
    opts.minSamplesPerWindow = 80;
    opts.PeakZ_warn = 1.0;
    opts.SNR_warn = 1.0;
    opts.nRMS_low = 0.10;
    opts.nRMS_high = 0.25;
    opts.dVar_low = 0.05;
    opts.dVar_high = 0.15;
    opts.SAVE_ANCHOR_FIGS = true;
    opts.SCALO_MIN_PERIOD_MIN = 60;
    opts.SCALO_MAX_PERIOD_MIN = 1590;
    opts.SCALO_YTICKS_H = 0:4:26;
    opts.SUBTRACT_FUNDAMENTAL_IN_SELECTIVE = cfg.hsub.subtractFundamentalInSelective;
end
