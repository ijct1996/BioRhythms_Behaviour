function run_extended_script5_transition_resync(extHandoffDir)
%RUN_EXTENDED_SCRIPT5_TRANSITION_RESYNC Extended Script 5 — transition resync + FDR.
%
%   Primary scientific engine for light-switch ultradian resynchronisation.
%   Co-primary bands: UR_1_3 and UR_3_6. Co-primary endpoints: DeltaR and
%   ridge-power change at DL/LD transitions.
%
%   Input: ExtendedHandoff/AcrossPhotoperiod_Input/ from Script 4.
%
%   run_extended_script5_transition_resync()
%   run_extended_script5_transition_resync(extHandoffDir)

    setup_extended_paths();
    cfg = extended_defaults();
    if nargin < 1
        extended_script5_run([], cfg);
    else
        extended_script5_run(extHandoffDir, cfg);
    end
end
