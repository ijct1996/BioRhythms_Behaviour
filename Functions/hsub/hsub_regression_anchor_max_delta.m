function bestDelta = hsub_regression_anchor_max_delta(t_hr, y, periods_h)
%HSUB_REGRESSION_ANCHOR_MAX_DELTA Max delta-R² for block-shuffle null.
%   Legacy port from harmonic_subtraction_v6.m

    t_hr = t_hr(:);
    y = y(:);
    [deltaCurve, ~, ~] = hsub_anchor_delta_curve(t_hr, y, periods_h);
    if all(~isfinite(deltaCurve))
        bestDelta = NaN;
    else
        bestDelta = max(deltaCurve, [], 'omitnan');
        if ~isfinite(bestDelta), bestDelta = NaN; end
    end
end
