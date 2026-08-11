function out = extended_script11_run(cohortRoot, cfg)
%EXTENDED_SCRIPT11_RUN Dominant validated UR period summary + supplemental violins.

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    cfg = extended_apply_plot_cfg(cfg);
    out = extended_script11_dominant_period_run(cohortRoot, cfg);
    fprintf('\nExtended Script 11 complete.\n');
end
