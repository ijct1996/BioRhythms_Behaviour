function run_extended_publication_profiles(varargin)
%RUN_EXTENDED_PUBLICATION_PROFILES Deprecated — use run_extended_script7_phase_profiles.
    warning('run_extended_publication_profiles:Deprecated', ...
        'Kent G is merged into Extended Script 7. Use run_extended_script7_phase_profiles.');
    setup_extended_paths();
    if nargin > 0
        run_extended_script7_phase_profiles(varargin{1});
    else
        run_extended_script7_phase_profiles();
    end
end
