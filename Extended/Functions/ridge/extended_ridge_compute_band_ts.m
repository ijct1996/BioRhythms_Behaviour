function bandTS = extended_ridge_compute_band_ts(periods_hours, logPow, bands, coi_hours, wt)
%EXTENDED_RIDGE_COMPUTE_BAND_TS Ridge-following band time series (Kent B logic).
    per = double(periods_hours(:));
    nT = size(logPow, 2);

    if isvector(bands) && numel(bands) == 2
        bands = reshape(bands, 1, 2);
    end
    nB = size(bands, 1);

    bandTS = struct();
    bandTS.BandPower = NaN(nB, nT);
    bandTS.RidgePeriod = NaN(nB, nT);
    bandTS.RidgePower = NaN(nB, nT);
    bandTS.RidgePhase = NaN(nB, nT);
    bandTS.ValidFlag = false(nB, nT);

    hasWT = (nargin >= 5) && ~isempty(wt) && isequal(size(wt), size(logPow));
    useCOI = ~isempty(coi_hours) && numel(coi_hours) >= nT;
    if useCOI
        coiH = double(coi_hours(:));
        coiH = coiH(1:nT);
    else
        coiH = [];
    end

    for b = 1:nB
        Pmin = bands(b, 1);
        Pmax = bands(b, 2);
        idx = find(per >= Pmin & per <= Pmax);
        if isempty(idx), continue; end

        sub = logPow(idx, :);
        validT = true(1, nT);
        if useCOI
            validT = isfinite(coiH).' & (coiH.' >= Pmax);
        end

        bp = mean(sub, 1, 'omitnan');
        [mx, k] = max(sub, [], 1, 'includenan');
        ok = validT & isfinite(bp) & isfinite(mx) & isfinite(k);
        bp(~ok) = NaN;
        mx(~ok) = NaN;
        k(~ok) = NaN;

        bandTS.BandPower(b, :) = bp;
        bandTS.RidgePower(b, :) = mx;
        perBand = per(idx);
        rp = NaN(1, nT);
        rp(ok) = perBand(k(ok));
        bandTS.RidgePeriod(b, :) = rp;

        if hasWT && any(ok)
            phaseRow = NaN(1, nT);
            okIdx = find(ok);
            kOK = k(ok);
            for jj = 1:numel(okIdx)
                try
                    phaseRow(okIdx(jj)) = angle(wt(idx(kOK(jj)), okIdx(jj)));
                catch
                    phaseRow(okIdx(jj)) = NaN;
                end
            end
            bandTS.RidgePhase(b, :) = phaseRow;
        end
        bandTS.ValidFlag(b, :) = ok;
    end
end
