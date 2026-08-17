function out = extended_script14_run(cohortRoot, cfg)
%EXTENDED_SCRIPT14_RUN Publication stitched scalograms (RAW + HSub residual, F|M).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    cfg = extended_apply_plot_cfg(cfg);
    out = extended_script14_scalogram_run(cohortRoot, cfg);
    fprintf('\nExtended Script 14 complete.\n');
end
