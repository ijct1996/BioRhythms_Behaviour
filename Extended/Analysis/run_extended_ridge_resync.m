function run_extended_ridge_resync(varargin)
%RUN_EXTENDED_RIDGE_RESYNC Deprecated — use run_extended_script5_transition_resync.
    warning('run_extended_ridge_resync:Deprecated', ...
        'Use run_extended_script5_transition_resync instead.');
    setup_extended_paths();
    if nargin > 0
        run_extended_script5_transition_resync(varargin{1});
    else
        run_extended_script5_transition_resync();
    end
end
