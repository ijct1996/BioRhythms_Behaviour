function run_extended_script9_supplementary_figures(cohortRoot)
%RUN_EXTENDED_SCRIPT9_SUPPLEMENTARY_FIGURES Extended Script 9 — supplementary figures.
%
%   Dedicated supplementary outputs (sex, primary-cluster 24h coherence + activity
%   L12–L22 grids, other-cluster activity). Main narrative figures remain Script 8.
%
%   run_extended_script9_supplementary_figures()
%   run_extended_script9_supplementary_figures(cohortRoot)

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 9 plot mode.\n\n' ...
        'Publication  - 600 DPI JPEG.\n' ...
        'Development  - 96 DPI PNG/JPEG for quick inspection.\n']), ...
        'Extended Script 9 Mode', ...
        'Publication', 'Development', 'Cancel', 'Publication');
    if isempty(choice) || strcmp(choice, 'Cancel')
        fprintf('Cancelled.\n');
        return;
    end
    if strcmp(choice, 'Publication')
        cfg.plotMode = 'publication';
    else
        cfg.plotMode = 'development';
    end
    cfg = extended_apply_plot_cfg(cfg);

    if nargin < 1
        extended_script9_run([], cfg);
    else
        extended_script9_run(cohortRoot, cfg);
    end
end
