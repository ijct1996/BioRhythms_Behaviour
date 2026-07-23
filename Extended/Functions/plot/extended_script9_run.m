function out = extended_script9_run(cohortRoot, cfg)
%EXTENDED_SCRIPT9_RUN Supplementary publication figures (Scripts 6–7 inputs).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    out = extended_script9_supplementary_figures_run(cohortRoot, cfg);
end
