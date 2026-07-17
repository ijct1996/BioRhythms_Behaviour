function run_extended_phase_events(varargin)
%RUN_EXTENDED_PHASE_EVENTS Deprecated — use run_extended_script7_phase_profiles.
    warning('run_extended_phase_events:Deprecated', ...
        'Kent F is merged into Extended Script 7. Use run_extended_script7_phase_profiles.');
    setup_extended_paths();
    if nargin > 0
        run_extended_script7_phase_profiles(varargin{1});
    else
        run_extended_script7_phase_profiles();
    end
end
