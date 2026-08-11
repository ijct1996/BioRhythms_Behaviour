function out = extended_script12_run(cohortRoot, cfg)
%EXTENDED_SCRIPT12_RUN Sex-stratified activity grids + post-LD amplitude (Script 12).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    cfg = extended_apply_plot_cfg(cfg);
    out = extended_script12_sex_stratified_run(cohortRoot, cfg);
    fprintf('\nExtended Script 12 complete.\n');
end
