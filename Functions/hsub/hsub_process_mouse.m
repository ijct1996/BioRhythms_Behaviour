function out = hsub_process_mouse(x0, timeMinutesAll, opts, LightDur_h)
%HSUB_PROCESS_MOUSE Harmonic subtraction for one mouse column.
%
%   Returns struct with residuals, reports, anchor info.

    out = struct();
    out.colName = '';
    out.anchorOK = false;
    out.P0_h = NaN;
    out.reports = struct();

    x0 = x0(:);
    missingOriginal = isnan(x0) | ~isfinite(x0);
    missingFrac = sum(missingOriginal) / numel(x0);
    maxGapSamples = max(0, round(opts.maxGapMinutes / opts.TsMinutes));
    doInterpolate = (missingFrac <= opts.missingFracThreshold) && (maxGapSamples > 0);

    xFill = x0;
    if doInterpolate
        xFill = fillmissing(x0, 'pchip', 'MaxGap', maxGapSamples);
    end

    base = movmedian(xFill, opts.baselineWinSamples, 'omitnan', 'Endpoints', 'shrink');
    yDetrFull = xFill - base;
    validIdx = find(~isnan(yDetrFull) & isfinite(timeMinutesAll));

    minPeriods = opts.minPeriodsToRun;
    keys = arrayfun(@(mp) sprintf('Min%d', mp), minPeriods, 'UniformOutput', false);

    for ki = 1:numel(keys)
        out.residuals.(keys{ki}).selective = x0;
        out.residuals.(keys{ki}).full = x0;
        out.residuals.(keys{ki}).removedSel = NaN(size(x0));
        out.residuals.(keys{ki}).removedFull = NaN(size(x0));
    end

    if numel(validIdx) < opts.minSamplesForAnalysis
        out.note = 'Too few valid samples';
        return;
    end

    t_min = timeMinutesAll(validIdx);
    t_hr = t_min / 60;
    durCol_h = (max(t_min, [], 'omitnan') - min(t_min, [], 'omitnan')) / 60;
    y = yDetrFull(validIdx);
    y = y - median(y, 'omitnan');

    periods_h = (opts.anchorBandHours(1):opts.anchorStepHours:opts.anchorBandHours(2))';
    [P0_h, maxDeltaR2, peakZ, A0, ampSNR, ~] = hsub_regression_anchor_scan(t_hr, y, periods_h);
    cyclesAtP0 = durCol_h / max(P0_h, eps);

    if numel(y) < 4 * opts.blockLenSamples
        blk = max(4, floor(numel(y) / 4));
    else
        blk = opts.blockLenSamples;
    end

    surrMax = NaN(opts.nSurrogates, 1);
    for s = 1:opts.nSurrogates
        yS = hsub_block_shuffle(y, blk);
        surrMax(s) = hsub_regression_anchor_max_delta(t_hr, yS, periods_h);
    end
    pBlock = (1 + sum(surrMax >= maxDeltaR2)) / (1 + opts.nSurrogates);

    anchorOK = isfinite(P0_h) && isfinite(maxDeltaR2) && isfinite(pBlock) && ...
        (pBlock < opts.alphaAnchor) && (cyclesAtP0 >= opts.minCyclesForAnchor) && ...
        (maxDeltaR2 >= opts.minDeltaR2);

    edgeOK = true;
    if opts.useEdgeGuard
        edgeOK = (P0_h > opts.anchorBandHours(1) + opts.edgeMarginHours) && ...
            (P0_h < opts.anchorBandHours(2) - opts.edgeMarginHours);
        anchorOK = anchorOK && edgeOK;
    end

    out.anchorOK = anchorOK;
    out.P0_h = P0_h;
    out.maxDeltaR2 = maxDeltaR2;
    out.peakZ = peakZ;
    out.pBlock = pBlock;
    out.cyclesAtP0 = cyclesAtP0;
    out.duration_h = durCol_h;
    out.TsMinutes = opts.TsMinutes;
    out.LightDur_h = LightDur_h;

    if ~anchorOK
        out.note = hsub_anchor_reject_note(P0_h, pBlock, cyclesAtP0, maxDeltaR2, ...
            opts, edgeOK, durCol_h);
        return;
    end

    winHours = min(max(72, round(durCol_h / 2)), 168);
    [winStarts_min, winEnds_min] = hsub_make_windows(min(t_min), max(t_min), winHours * 60, opts.stepHours * 60);
    nW = numel(winStarts_min);
    P0_min = P0_h * 60;

    A1_w = NaN(nW, 1);
    Phi1_w = NaN(nW, 1);
    for w = 1:nW
        wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
        if sum(wMask) < opts.minSamplesPerWindow, continue; end
        [A1, phi1, ~] = hsub_fit_sinusoid(t_min(wMask), y(wMask), P0_min);
        A1_w(w) = A1;
        Phi1_w(w) = phi1;
    end

    okFund = isfinite(A1_w) & isfinite(Phi1_w);
    nFund = sum(okFund);
    if nFund >= 4
        regime = 'long';
        minWindowsForDecision = 4;
    elseif nFund >= 2
        regime = 'short';
        minWindowsForDecision = 2;
    else
        regime = 'insufficient';
        minWindowsForDecision = inf;
    end

    crossQC = struct();

    for mp = minPeriods
        key = sprintf('Min%d', mp);
        K = max(2, min(24, floor(P0_min / mp)));
        kList = 2:K;

        Ak_w = NaN(nW, numel(kList));
        Phik_w = NaN(nW, numel(kList));
        ratio_w = NaN(nW, numel(kList));

        for w = 1:nW
            wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
            if sum(wMask) < opts.minSamplesPerWindow || ~isfinite(A1_w(w)), continue; end
            tmw = t_min(wMask);
            yw = y(wMask);
            for j = 1:numel(kList)
                k = kList(j);
                targetP = P0_min / k;
                scanP = linspace(targetP * (1 - opts.peakSearchFrac), targetP * (1 + opts.peakSearchFrac), 19);
                bestR2 = -Inf;
                bestA = NaN; bestPhi = NaN; bestP = NaN;
                for sp = 1:numel(scanP)
                    [Akk, phikk, R2kk] = hsub_fit_sinusoid(tmw, yw, scanP(sp));
                    if isfinite(R2kk) && R2kk > bestR2
                        bestR2 = R2kk; bestA = Akk; bestPhi = phikk; bestP = scanP(sp);
                    end
                end
                if ~isfinite(bestR2), continue; end
                Ak_w(w, j) = bestA;
                Phik_w(w, j) = bestPhi;
                ratio_w(w, j) = bestP / targetP;
            end
        end

        harmonicLikely = false(1, numel(kList));
        for j = 1:numel(kList)
            okW = isfinite(A1_w) & isfinite(Ak_w(:, j)) & isfinite(Phi1_w) & isfinite(Phik_w(:, j));
            if sum(okW) < minWindowsForDecision, continue; end
            ratios = ratio_w(okW, j);
            ratioPass = mad(ratios - 1, 1) <= opts.ratioTol;
            del = hsub_wrap_pi(Phik_w(okW, j) - kList(j) * Phi1_w(okW));
            Rplv = abs(mean(exp(1i * del)));
            Rpass = isfinite(Rplv) && (Rplv >= opts.Rcrit);
            if strcmp(regime, 'long')
                rho = corr(A1_w(okW), Ak_w(okW, j), 'Type', 'Spearman', 'Rows', 'complete');
                harmonicLikely(j) = (ratioPass + (isfinite(rho) && abs(rho) >= opts.rhoCrit) + Rpass) >= 2;
            else
                A1_med = median(A1_w(okW), 'omitnan');
                Ak_med = median(Ak_w(okW, j), 'omitnan');
                ampPass = isfinite(A1_med) && A1_med > 0 && isfinite(Ak_med) && ...
                    (Ak_med >= opts.ampFracCritShort * A1_med);
                harmonicLikely(j) = ratioPass && Rpass && ampPass;
            end
        end

        kLikely = kList(harmonicLikely);
        if opts.SUBTRACT_FUNDAMENTAL_IN_SELECTIVE
            kSubSel = unique([1, kLikely], 'stable');
        else
            kSubSel = kLikely;
        end

        if isempty(kSubSel)
            residualSel = y;
            removedSel = zeros(size(y));
            varExplSel = 0;
        else
            [X0, Xh] = hsub_build_design(t_hr, t_min, P0_min, kSubSel);
            beta = [X0, Xh] \ y;
            yhatDrift = X0 * beta(1:size(X0, 2));
            yhatTotal = [X0, Xh] * beta;
            removedSel = yhatTotal - yhatDrift;
            residualSel = y - removedSel;
            varExplSel = hsub_var_explained(y, residualSel);
        end

        kSubFull = 1:K;
        [X0F, XhF] = hsub_build_design(t_hr, t_min, P0_min, kSubFull);
        betaF = [X0F, XhF] \ y;
        yhatDriftF = X0F * betaF(1:size(X0F, 2));
        yhatTotalF = [X0F, XhF] * betaF;
        removedFull = yhatTotalF - yhatDriftF;
        residualFull = y - removedFull;
        varExplFull = hsub_var_explained(y, residualFull);

        [rSel, remSel] = hsub_map_full(x0, validIdx, residualSel, removedSel);
        [rFull, remFull] = hsub_map_full(x0, validIdx, residualFull, removedFull);

        out.residuals.(key).selective = rSel;
        out.residuals.(key).full = rFull;
        out.residuals.(key).removedSel = remSel;
        out.residuals.(key).removedFull = remFull;
        out.residuals.(key).varExplSel = varExplSel;
        out.residuals.(key).varExplFull = varExplFull;
        out.residuals.(key).K = K;

        rmsDiff = sqrt(mean((residualSel - residualFull).^2, 'omitnan'));
        nRMS = rmsDiff / (sqrt(mean(y.^2, 'omitnan')) + eps);
        crossQC.(key).nRMS = nRMS;
        crossQC.(key).dVar = varExplFull - varExplSel;
    end

    out.crossQC = crossQC;
    out.note = sprintf('OK (%s regime)', regime);
end

function yS = hsub_block_shuffle(y, blockLen)
    y = y(:); n = numel(y);
    if n < 2 * blockLen, yS = y(randperm(n)); return; end
    nB = floor(n / blockLen); remN = n - nB * blockLen;
    blocks = cell(nB + (remN > 0), 1);
    for b = 1:nB
        i1 = (b-1) * blockLen + 1; i2 = b * blockLen;
        blocks{b} = y(i1:i2);
    end
    if remN > 0, blocks{end} = y(nB * blockLen + 1:end); end
    yS = vertcat(blocks{randperm(numel(blocks))});
end

function [A, phi, R2] = hsub_fit_sinusoid(t_min, y, period_min)
    t = t_min(:); y = y(:);
    if numel(t) < 10 || all(~isfinite(y)), A = NaN; phi = NaN; R2 = NaN; return; end
    w = 2 * pi / period_min;
    tZ = hsub_zscore(t);
    X = [cos(w*t), sin(w*t), ones(size(t)), tZ];
    b = X \ y; yhat = X * b;
    A = hypot(b(1), b(2)); phi = atan2(-b(2), b(1));
    ssRes = sum((y - yhat).^2, 'omitnan');
    ssTot = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
    R2 = 1 - ssRes / (ssTot + eps);
end

function [X0, Xh] = hsub_build_design(t_hr, t_min, P0_min, kList)
    tz = hsub_zscore(t_hr(:));
    X0 = [ones(size(tz)), tz, tz.^2];
    t = t_min(:); nK = numel(kList);
    Xh = zeros(numel(t), 2 * nK);
    for ii = 1:nK
        w = 2 * pi * kList(ii) / P0_min;
        col = 2 * (ii - 1) + 1;
        Xh(:, col) = cos(w * t);
        Xh(:, col + 1) = sin(w * t);
    end
end

function z = hsub_zscore(x)
    x = x(:); sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0, z = x - mean(x, 'omitnan'); else, z = (x - mean(x, 'omitnan')) / sd; end
end

function x = hsub_wrap_pi(x), x = mod(x + pi, 2*pi) - pi; end

function [starts, ends] = hsub_make_windows(tMin, tMax, winLenMin, stepMin)
    starts = (tMin:stepMin:(tMax - winLenMin))';
    ends = starts + winLenMin;
    if isempty(starts), starts = tMin; ends = tMax; end
end

function [rFull, remFull] = hsub_map_full(x0, validIdx, residual, removed)
    rFull = x0; remFull = NaN(size(x0));
    rFull(validIdx) = residual; remFull(validIdx) = removed;
    rFull(isnan(x0)) = NaN; remFull(isnan(x0)) = NaN;
end

function v = hsub_var_explained(y, residual)
    vy = var(y, 0, 'omitnan');
    if ~isfinite(vy) || vy <= 0, v = NaN; else, v = 1 - var(residual, 0, 'omitnan') / (vy + eps); end
end

function note = hsub_anchor_reject_note(P0_h, pBlock, cyclesAtP0, maxDeltaR2, opts, edgeOK, durCol_h)
    reasons = {};
    if ~isfinite(P0_h) || ~isfinite(maxDeltaR2) || ~isfinite(pBlock)
        reasons{end+1} = 'Non-finite anchor statistics';
    end
    if isfinite(pBlock) && pBlock >= opts.alphaAnchor
        reasons{end+1} = sprintf('pBlock=%.3g (>=%.2f)', pBlock, opts.alphaAnchor);
    end
    if isfinite(cyclesAtP0) && cyclesAtP0 < opts.minCyclesForAnchor
        reasons{end+1} = sprintf('cycles=%.2f at P0 (need >=%.1f; duration=%.1f h)', ...
            cyclesAtP0, opts.minCyclesForAnchor, durCol_h);
    end
    if isfinite(maxDeltaR2) && maxDeltaR2 < opts.minDeltaR2
        reasons{end+1} = sprintf('maxDeltaR2=%.4g (<%.4g)', maxDeltaR2, opts.minDeltaR2);
    end
    if ~edgeOK
        reasons{end+1} = sprintf('Period=%.3f h too close to band edge', P0_h);
    end
    if isempty(reasons)
        note = 'Anchor rejected';
    else
        note = strjoin(reasons, '; ');
    end
end
