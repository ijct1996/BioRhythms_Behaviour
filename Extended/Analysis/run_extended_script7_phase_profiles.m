function run_extended_script7_phase_profiles(extHandoffDir)
%RUN_EXTENDED_SCRIPT7_PHASE_PROFILES Extended Script 7 — phase events + profiles.
%
%   Local Extended implementation with development/publication mode chooser.
%
%   run_extended_script7_phase_profiles()
%   run_extended_script7_phase_profiles(extHandoffDir)

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 7 plot mode.\n\n' ...
        'Publication  - high DPI outputs for manuscript-ready figures.\n' ...
        'Development  - faster/lighter outputs for inspection.\n']), ...
        'Extended Script 7 Mode', ...
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
        extended_script7_run([], cfg);
    else
        extended_script7_run(extHandoffDir, cfg);
    end
end
