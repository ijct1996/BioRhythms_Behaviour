function run_extended_script12_sex_stratified_profiles(cohortRoot)
%RUN_EXTENDED_SCRIPT12_SEX_STRATIFIED_PROFILES Extended Script 12 — activity + pre/post-LD amplitude.
%
%   Outputs for UR_1_3 C01 | UR_3_6 C01 | UR_3_6 C02:
%     - Female|Male side-by-side 24h activity grids (L12–L22)
%     - Pre– and post–lights-off first-peak amplitude vs photoperiod (F/M + BH-FDR stars)
%
%   No coherence / ridge-power panels. Read-only Script 7 profiles.
%
%   run_extended_script12_sex_stratified_profiles()
%   run_extended_script12_sex_stratified_profiles(cohortRoot)

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 12 plot mode.\n\n' ...
        'Publication  - 600 DPI JPEG (+ TIFF when enabled).\n' ...
        'Development  - lighter PNG/JPEG for inspection.\n']), ...
        'Extended Script 12 Mode', ...
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
        extended_script12_run([], cfg);
    else
        extended_script12_run(cohortRoot, cfg);
    end
end
