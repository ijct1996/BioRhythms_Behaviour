function run_extended_script14_scalograms(cohortRoot)
%RUN_EXTENDED_SCRIPT14_SCALOGRAMS Extended Script 14 — publication stitched scalograms.
%
%   Recomputes sex-split stitched group-mean CWT scalograms from Handoff/ +
%   HSub TimeSeries (SEL_P360 residual). Four outputs: RAW F|M, HSub residual F|M.
%   Sex-pooled color limits within RAW ([0, 14]) and HSub ([0, 10]).
%
%   run_extended_script14_scalograms()
%   run_extended_script14_scalograms(cohortRoot)
%
%   Non-interactive:
%     cfg = extended_defaults(); cfg.plotMode = 'publication';
%     cfg = extended_apply_plot_cfg(cfg);
%     extended_script14_run(cohortRoot, cfg)

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 14 plot mode.\n\n' ...
        'Publication  - 600 DPI JPEG scalograms (Arial).\n' ...
        'Development  - 96 DPI PNG for quick inspection.\n']), ...
        'Extended Script 14 Mode', ...
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
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);
    fprintf('Script 14 mode locked: %s (dpi=%g, ext=%s)\n', ...
        cfg.plotMode, cfg.plot.saveDpi, char(string(cfg.plot.figExt)));

    if nargin < 1
        extended_script14_run([], cfg);
    else
        extended_script14_run(cohortRoot, cfg);
    end
end
