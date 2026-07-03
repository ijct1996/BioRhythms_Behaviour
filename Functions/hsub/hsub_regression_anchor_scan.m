function [P0_h, maxDeltaR2, peakZ, amp, ampSNR, deltaCurve] = hsub_regression_anchor_scan(t_hr, y, periods_h)
%HSUB_REGRESSION_ANCHOR_SCAN Legacy regression anchor scan (22–28 h band).
%
%   Ported from harmonic_subtraction_v6.m — period in HOURS, time t_hr in HOURS.

    t_hr = t_hr(:);
    y = y(:);

    [deltaCurve, ampCurve, y0] = hsub_anchor_delta_curve(t_hr, y, periods_h);

    if all(~isfinite(deltaCurve))
        P0_h = NaN; maxDeltaR2 = NaN; peakZ = NaN; amp = NaN; ampSNR = NaN;
        return;
    end

    [maxDeltaR2, iBest] = max(deltaCurve, [], 'omitnan');
    if ~isfinite(maxDeltaR2) || isempty(iBest) || iBest < 1 || iBest > numel(periods_h)
        P0_h = NaN; amp = NaN;
    else
        P0_h = periods_h(iBest);
        amp = ampCurve(iBest);
    end

    medD = median(deltaCurve, 'omitnan');
    madD = mad(deltaCurve, 1);
    if ~isfinite(madD) || madD == 0
        peakZ = Inf;
    else
        peakZ = (maxDeltaR2 - medD) / (madD + eps);
    end

    r0 = y - y0;
    robustStd = 1.4826 * mad(r0, 1);
    if ~isfinite(robustStd) || robustStd <= 0 || ~isfinite(amp)
        ampSNR = 0;
    else
        ampSNR = amp / robustStd;
    end
end
