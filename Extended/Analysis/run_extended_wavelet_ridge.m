function run_extended_wavelet_ridge(varargin)
%RUN_EXTENDED_WAVELET_RIDGE Deprecated — use run_extended_script4_ridge_validation.
    warning('run_extended_wavelet_ridge:Deprecated', ...
        'Kent B ridge handoff is in Extended Script 4. Use run_extended_script4_ridge_validation.');
    setup_extended_paths();
    if nargin > 0
        run_extended_script4_ridge_validation(varargin{1});
    else
        run_extended_script4_ridge_validation();
    end
end
