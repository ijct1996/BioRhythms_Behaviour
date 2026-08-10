function out = extended_script10_run(cohortRoot, cfg)
%EXTENDED_SCRIPT10_RUN Manuscript figure assembler (Scripts 5–9 inputs).
%
%   Non-interactive entry: pass cfg.plotMode ('publication'|'development').
%   Prefer run_extended_script10_manuscript_figures for interactive mode dialog.

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    % Lock mode before defaults/theme so folder + dpi cannot silently fall back
    if ~isfield(cfg, 'plotMode') || isempty(cfg.plotMode)
        cfg.plotMode = 'development';
    end
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);

    out = extended_script10_manuscript_figures_run(cohortRoot, cfg);
end
