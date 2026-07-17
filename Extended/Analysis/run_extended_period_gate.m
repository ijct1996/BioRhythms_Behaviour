function run_extended_period_gate(varargin)
%RUN_EXTENDED_PERIOD_GATE Deprecated — use run_extended_script4_ridge_validation.
    warning('run_extended_period_gate:Deprecated', ...
        'Use run_extended_script4_ridge_validation instead.');
    setup_extended_paths();
    if nargin > 0
        run_extended_script4_ridge_validation(varargin{1});
    else
        run_extended_script4_ridge_validation();
    end
end
