function signal = extended_wavelet_prepare_vector(x)
%EXTENDED_WAVELET_PREPARE_VECTOR Numeric vector with non-finite → 0.
    signal = x(:);
    if ~isnumeric(signal)
        signal = str2double(string(signal));
    end
    signal(~isfinite(signal)) = 0;
end
