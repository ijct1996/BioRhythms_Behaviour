function signal = wavelet_prepare_signal(tbl, colIdx)
%WAVELET_PREPARE_SIGNAL Numeric column vector with non-finite → 0.
    signal = tbl{:, colIdx};
    if ~isnumeric(signal)
        signal = str2double(string(signal));
    end
    signal = signal(:);
    signal(~isfinite(signal)) = 0;
end
