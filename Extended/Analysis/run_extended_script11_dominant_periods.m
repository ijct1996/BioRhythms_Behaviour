function run_extended_script11_dominant_periods(cohortRoot)
%RUN_EXTENDED_SCRIPT11_DOMINANT_PERIODS Extended Script 11 — dominant UR period summary + supp violins.
%
%   Read-only post-hoc: requires Script 4 (WP_TS) + Script 7 (cluster tables).
%   Produces manuscript median-period tables and a 3-panel supplemental figure
%   (UR_1_3 C01 | UR_3_6 C01 | UR_3_6 C02).
%
%   run_extended_script11_dominant_periods()
%   run_extended_script11_dominant_periods(cohortRoot)

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 11 plot mode.\n\n' ...
        'Publication  - 600 DPI JPEG (+ TIFF when enabled).\n' ...
        'Development  - lighter PNG/JPEG for inspection.\n']), ...
        'Extended Script 11 Mode', ...
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
        extended_script11_run([], cfg);
    else
        extended_script11_run(cohortRoot, cfg);
    end
end
