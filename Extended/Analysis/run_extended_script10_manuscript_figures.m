function run_extended_script10_manuscript_figures(cohortRoot)
%RUN_EXTENDED_SCRIPT10_MANUSCRIPT_FIGURES Extended Script 10 — locked manuscript figures.
%
%   Assembles Scripts 5–9 outputs into the locked main figure set (Figs 1–7).
%   Figs 5 & 7 main composites = 24h ZT coherence (A) only; Script 5 B/C
%   (ridge power / DeltaR) export under Standalone/Supp_Transitions/.
%   UR_3_6 C02 activity+coherence twins under Standalone/Supp/ (cfg lock).
%   Builds Package_Manuscript/ (Figs 3–7 + Figures/Supp + tables; excludes Fig 1 & 2).
%
%   run_extended_script10_manuscript_figures()
%   run_extended_script10_manuscript_figures(cohortRoot)
%
%   Non-interactive (batch-safe):
%     cfg = extended_defaults(); cfg.plotMode = 'development';
%     cfg = extended_apply_plot_cfg(cfg);
%     extended_script10_run(cohortRoot, cfg);

    setup_extended_paths();
    cfg = extended_defaults();
    choice = questdlg(sprintf([ ...
        'Select Script 10 plot mode.\n\n' ...
        'Publication  - 600 DPI JPEG + PDF vectors for line plots.\n' ...
        'Development  - 96 DPI PNG for quick inspection.\n']), ...
        'Extended Script 10 Mode', ...
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
    % Force char + re-apply so dpi/ext/folder cannot stick on defaults('development')
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);
    fprintf('Script 10 mode locked: %s (dpi=%g, ext=%s)\n', ...
        cfg.plotMode, cfg.plot.saveDpi, char(string(cfg.plot.figExt)));

    if nargin < 1
        extended_script10_run([], cfg);
    else
        extended_script10_run(cohortRoot, cfg);
    end
end
