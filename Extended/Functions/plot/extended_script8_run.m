function out = extended_script8_run(cohortRoot, cfg)
%EXTENDED_SCRIPT8_RUN Publication figure composites from Scripts 1-7 outputs.

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    out = extended_script8_publication_figures_run(cohortRoot, cfg);
end
