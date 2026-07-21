function run_extended_script8_publication_figures(cohortRoot)
%RUN_EXTENDED_SCRIPT8_PUBLICATION_FIGURES Extended Script 8 — publication composites.
%
%   Select cohort results folder (e.g. C57_LP). Reads Scripts 1-7 outputs in place.
%
%   run_extended_script8_publication_figures()
%   run_extended_script8_publication_figures(cohortRoot)

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 8 plot mode.\n\n' ...
        'Publication  - 600 DPI JPEG + PDF vectors for line plots.\n' ...
        'Development  - 96 DPI PNG for quick inspection.\n']), ...
        'Extended Script 8 Mode', ...
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
        extended_script8_run([], cfg);
    else
        extended_script8_run(cohortRoot, cfg);
    end
end
