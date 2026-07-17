function run_extended_script6_across_lme(extHandoffDir)
%RUN_EXTENDED_SCRIPT6_ACROSS_LME Extended Script 6 — LME / FDR across photoperiod.
%
%   Complements Core Script 3 descriptive outputs with inferential LME on
%   CarryForward-validated raw UR band power (Kent E).
%
%   run_extended_script6_across_lme()
%   run_extended_script6_across_lme(extHandoffDir)

    setup_extended_paths();
    cfg = extended_defaults();
    if nargin < 1
        extended_script6_run([], cfg);
    else
        extended_script6_run(extHandoffDir, cfg);
    end
end
