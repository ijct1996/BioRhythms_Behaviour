function signalsMatrix = wavelet_stack_signals(tbl, colIdx)
%WAVELET_STACK_SIGNALS Matrix of activity columns for group averaging.
    signalsMatrix = [];
    for k = 1:numel(colIdx)
        signalsMatrix(:, k) = wavelet_prepare_signal(tbl, colIdx(k)); %#ok<AGROW>
    end
end
