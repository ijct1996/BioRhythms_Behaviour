function out = extended_script6_run(extHandoffDir, cfg)
%EXTENDED_SCRIPT6_RUN Across-photoperiod LME / FDR on validated raw UR (local).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    out = extended_script6_across_lme_run(extHandoffDir, cfg);
end
