function [deltaCurve, ampCurve, y0] = hsub_anchor_delta_curve(t_hr, y, periods_h)
%HSUB_ANCHOR_DELTA_CURVE Drift + sinusoid scan across candidate periods (hours).
%   Legacy port from harmonic_subtraction_v6.m

    t_hr = t_hr(:);
    y = y(:);

    tz = hsub_zscore_local(t_hr);
    X0 = [ones(size(t_hr)), tz, tz.^2];
    b0 = X0 \ y;
    y0 = X0 * b0;

    ssTot = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
    ssRes0 = sum((y - y0).^2, 'omitnan');
    R2_0 = 1 - (ssRes0 / (ssTot + eps));

    deltaCurve = NaN(size(periods_h));
    ampCurve = NaN(size(periods_h));

    for i = 1:numel(periods_h)
        P = periods_h(i);
        w = 2 * pi / P;
        X = [X0, cos(w * t_hr), sin(w * t_hr)];
        b = X \ y;
        yhat = X * b;

        ssRes = sum((y - yhat).^2, 'omitnan');
        R2 = 1 - (ssRes / (ssTot + eps));
        deltaCurve(i) = R2 - R2_0;

        a = b(end - 1);
        s = b(end);
        ampCurve(i) = hypot(a, s);
    end
end

function z = hsub_zscore_local(x)
    x = x(:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        z = x - mu;
    else
        z = (x - mu) / sd;
    end
end
