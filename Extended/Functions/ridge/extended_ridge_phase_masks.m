function [isLight, lightStateStr, phaseMasks] = extended_ridge_phase_masks(time_hr, lightDurVec)
%EXTENDED_RIDGE_PHASE_MASKS Light/dark masks from light-duration column (Kent rule).
    t = double(time_hr(:));
    L = double(lightDurVec(:));
    n = min(numel(t), numel(L));
    t = t(1:n);
    L = L(1:n);

    ZT = mod(t, 24);
    isLight = false(n, 1);
    ok = isfinite(ZT) & isfinite(L) & L >= 0 & L <= 24;
    isLight(ok) = (ZT(ok) < L(ok));

    ls = repmat("Dark", n, 1);
    ls(isLight) = "Light";
    lightStateStr = cellstr(ls);

    phaseMasks = struct();
    phaseMasks.All = true(n, 1);
    phaseMasks.Light = isLight(:);
    phaseMasks.Dark = ~phaseMasks.Light;
end
