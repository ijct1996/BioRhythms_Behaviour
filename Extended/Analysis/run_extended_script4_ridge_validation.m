function run_extended_script4_ridge_validation(coreHandoffDir)
%RUN_EXTENDED_SCRIPT4_RIDGE_VALIDATION Extended Script 4 — ridge handoff + gate.
%
%   Reads Core Script 1–3 Handoff/ (CoreSummary__*.mat). Adds ridge tracking
%   on RAW + HSub residuals only; runs CarryForward validation (±15% SEL_P360).
%   Does NOT duplicate Core HSub or scalogram wavelet outputs.
%
%   run_extended_script4_ridge_validation()
%   run_extended_script4_ridge_validation(coreHandoffDir)
%
%   See also: extended_script4_run, setup_extended_paths

    setup_extended_paths();
    cfg = extended_defaults();
    if nargin < 1
        extended_script4_run([], cfg);
    else
        extended_script4_run(coreHandoffDir, cfg);
    end
end
