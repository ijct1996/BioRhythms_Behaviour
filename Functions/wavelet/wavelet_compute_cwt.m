function [wt, periods_hours, coi_hours] = wavelet_compute_cwt(signal, FB)
%WAVELET_COMPUTE_CWT Continuous wavelet transform.
    signal = signal(:);
    signal(~isfinite(signal)) = 0;
    [wt, periods, coi] = cwt(signal, 'FilterBank', FB);
    periods_hours = hours(periods);
    coi_hours = hours(coi);
end
