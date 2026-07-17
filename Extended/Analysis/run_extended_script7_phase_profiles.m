function run_extended_script7_phase_profiles(extHandoffDir)
%RUN_EXTENDED_SCRIPT7_PHASE_PROFILES Extended Script 7 — phase events + profiles.
%
%   Merged Kent F (phase-event histograms) and G (publication profiles).
%   Skipped in development when cfg.plot.generateLegacyFigures=false.
%
%   run_extended_script7_phase_profiles()
%   run_extended_script7_phase_profiles(extHandoffDir)

    setup_extended_paths();
    cfg = extended_defaults();
    if nargin < 1
        extended_script7_run([], cfg);
    else
        extended_script7_run(extHandoffDir, cfg);
    end
end
