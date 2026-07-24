% -------------------------------------------------------------------------
% LB_01_HarmonicSubtraction.m
% -------------------------------------------------------------------------
% Circadian harmonic subtraction for locomotor activity time series.
%
% Developed for MATLAB R2025b. All helper functions used by this analysis
% sit at the bottom of this file.
%
% What this code does (short version)
%   1. Reads activity Excel files (Time, Light duration, mouse columns).
%   2. For each mouse, finds a circadian anchor period in 22-28 h.
%   3. Tests that anchor against block-shuffle surrogates.
%   4. If the anchor passes, subtracts the fundamental and likely harmonics
%      (selective and full-ladder arms; Min360 and Min60 period floors).
%   5. Writes residuals, removed components, and an Anchor/Recommendation report.
%   6. Writes SEL_P360 Removed|Residual scalograms (Figures_Scalograms/;
%      AnchorOK mice only). Wavelet Toolbox required for that step.
%
% Default recommended residual when the anchor is OK: SEL_P360
%   (Selective subtraction, minimum period floor 360 min).
%
% Sampling assumed: 10 min (0.167 h), unless Time (hr) implies otherwise.
% -------------------------------------------------------------------------

clearvars; close all; clc;
rng(1);   % fixed seed for surrogate reproducibility when re-running

fprintf('Harmonic subtraction\n\n');

%% ========================= USER CHOICES (DIALOGS) ========================
% Everything interactive happens here. After this block the run is unsupervised.

% Allows user to choose either single or multiple files - useful for dev
% (single files) and streamlining (multiple files)
runMode = questdlg( ...
    'Choose input mode:', ...
    'Harmonic subtraction', ...
    'Single file', 'Multiple files', 'Cancel', 'Single file');
if isempty(runMode) || strcmpi(runMode, 'Cancel')
    fprintf('Cancelled.\n');
    return;
end

[fileList, inPath] = lb_select_input_files(runMode);
if isempty(fileList)
    fprintf('No input file(s) selected. Exiting.\n');
    return;
end

% Prompts user to choose the output filepath
outRoot = uigetdir(inPath, 'Select output root (HarmonicSubtraction will be created here if needed)');
if isequal(outRoot, 0)
    fprintf('No output folder selected. Exiting.\n');
    return;
end

% Optional: nest under a standard folder name
if ~contains(lower(outRoot), 'harmonicsubtraction')
    outRoot = fullfile(outRoot, 'HarmonicSubtraction');
end
if ~exist(outRoot, 'dir'), mkdir(outRoot); end

plotMode = questdlg('Figure export mode (scalograms):', 'Plot mode', ...
    'development', 'publication', 'development');
if isempty(plotMode), plotMode = 'development'; end
% development = quick PNG (~96 DPI); publication = JPEG (~600 DPI)

saveIndivAns = questdlg( ...
    'Also save individual Removed|Residual scalograms (AnchorOK mice)?', ...
    'HSub scalograms', 'No', 'Yes', 'No');
% Group averages are always written; individuals are optional QC figures
saveIndividualScalograms = strcmpi(saveIndivAns, 'Yes');

opts = lb_hsub_default_opts();
fprintf('Output root: %s\n', outRoot);
fprintf('Files: %d | Anchor band: %.0f-%.0f h | Surrogates: %d | alpha: %.2f\n\n', ...
    numel(fileList), opts.anchorBandHours(1), opts.anchorBandHours(2), ...
    opts.nSurrogates, opts.alphaAnchor);

%% ========================= BATCH OVER INPUT FILES ========================
% One unsupervised pass per selected Excel file (errors are caught per file).
% Condition groups are defined on the first file and reused thereafter
% (needed for Figures_Scalograms/ group averages).

groupsTemplate = [];
for iF = 1:numel(fileList)
    inputFile = fileList{iF};
    [~, fileStem, ~] = fileparts(inputFile);
    fprintf('============================================================\n');
    fprintf('File %d/%d: %s\n', iF, numel(fileList), fileStem);
    fprintf('============================================================\n');

    try
        groupsTemplate = lb_hsub_run_one_file(inputFile, outRoot, opts, ...
            groupsTemplate, plotMode, saveIndividualScalograms);
    catch ME
        fprintf(2, 'FAILED: %s\n%s\n', fileStem, ME.message);
    end
end

fprintf('\nHarmonic subtraction batch finished.\n');
fprintf('Per-file folders are under:\n  %s\n', outRoot);
fprintf('Look for Reports/, TimeSeries/, and Figures_Scalograms/.\n');

%% ========================================================================
% LOCAL FUNCTIONS — harmonic subtraction
% (Reading order: run_one_file → process_mouse → math helpers)
% ========================================================================

function opts = lb_hsub_default_opts()
% Analysis defaults for harmonic subtraction.
% Change these only if you are deliberately changing the analysis protocol.
    opts.TsMinutes = 10;                 % expected sampling step (minutes)
    opts.minPeriodsMinutes = [360, 60];  % two arms: 6 h floor and 1 h floor
    opts.anchorBandHours = [22, 28];     % circadian search band for P0
    opts.anchorStepHours = 0.01;         % fine grid over that band
    opts.blockLenHours = 8;              % surrogate block length
    opts.nSurrogates = 300;
    opts.alphaAnchor = 0.05;             % reject anchor if pBlock >= this
    opts.missingFracThreshold = 0.05;    % only fill gaps if missingness is low
    opts.maxGapMinutes = 30;
    opts.minCyclesForAnchor = 3.5;       % need enough cycles at P0
    opts.minDeltaR2 = 0.005;             % tiny improvements are not accepted
    opts.useEdgeGuard = true;            % reject peaks stuck on 22 or 28 h
    opts.edgeMarginHours = 0.10;
    opts.subtractFundamentalInSelective = true;
    opts.peakSearchFrac = 0.06;          % ±6% search around P0/k for each harmonic
    opts.ratioTol = 0.08;                % period-ratio consistency for "likely" harmonic
    opts.Rcrit = 0.55;                   % phase-locking (PLV) threshold
    opts.rhoCrit = 0.45;                 % amplitude correlation (long recordings)
    opts.ampFracCritShort = 0.04;        % relative amplitude rule (short recordings)
    opts.stepHours = 24;                 % sliding-window step
    opts.minSamplesForAnalysis = 250;
    opts.minSamplesPerWindow = 80;
    % Sensitivity options below 
    opts.nRMS_low = 0.10;                % selective vs full residual QC bands
    opts.nRMS_high = 0.25;
    opts.dVar_low = 0.05;
    opts.dVar_high = 0.15;
    opts.defaultResidual = 'SEL_P360';
    % Scalogram export (SEL_P360 Removed|Residual after the residuals are written).
    % Morlet CWT settings match the wavelet analysis code so figures are comparable.
    opts.scalogramArm = 'Min360';
    opts.scalogramLabel = 'SEL_P360';
    opts.waveletName = 'amor';
    opts.periodLimitsMinutes = [60, 1590];
    opts.scalogramYTicksHours = 0:4:26;
    opts.colormap = 'jet';
    opts.plot.development = struct('dpi', 96, 'format', 'png');
    opts.plot.publication = struct('dpi', 600, 'format', 'jpeg');
    opts.indivScalogram = struct('dpi', 150, 'format', 'png');
end

function [fileList, inPath] = lb_select_input_files(runMode)
% Pop-up dialog(s) to choose one Excel/CSV file, or several from a folder.
% Returns fileList = cell array of full paths, and inPath = folder used.
    fileList = {};
    inPath = pwd;   % start looking from the current MATLAB folder
    if strcmpi(runMode, 'Single file')
        [f, p] = uigetfile({'*.xlsx;*.xls;*.csv', 'Excel / CSV'}, 'Select activity Excel file');
        if isequal(f, 0), return; end   % user pressed Cancel
        inPath = p;
        fileList = {fullfile(p, f)};    % one path inside a cell
    else
        inPath = uigetdir(pwd, 'Select folder containing activity Excel files');
        if isequal(inPath, 0), inPath = pwd; return; end
        % List .xlsx and .csv in that folder
        d = [dir(fullfile(inPath, '*.xlsx')); dir(fullfile(inPath, '*.csv'))];
        d = d(~[d.isdir]);
        if isempty(d)
            fprintf('No .xlsx/.csv files found in that folder.\n');
            return;
        end
        names = unique({d.name}, 'stable');
        [sel, ok] = listdlg('PromptString', 'Select files to process:', ...
            'SelectionMode', 'multiple', 'ListString', names, 'InitialValue', 1:numel(names));
        if ~ok || isempty(sel), return; end
        fileList = fullfile(inPath, names(sel));
        if ischar(fileList), fileList = {fileList}; end
    end
end

function groupsTemplate = lb_hsub_run_one_file(inputFile, outputRoot, opts, ...
        groupsTemplate, plotMode, saveIndividualScalograms)
% Read one Excel file, map columns, run every mouse, write residuals + report
% + SEL_P360 Removed|Residual scalograms under Figures_Scalograms/.

    if nargin < 4, groupsTemplate = []; end
    if nargin < 5 || isempty(plotMode), plotMode = 'development'; end
    if nargin < 6, saveIndividualScalograms = false; end

    % --- Import and column mapping ---
    [tbl, varNames] = lb_read_excel(inputFile);
    [timeIdx, lightIdx, mouseIdx] = lb_map_columns(varNames);
    if isempty(timeIdx) || isempty(lightIdx) || isempty(mouseIdx)
        error('Could not resolve Time / Light duration / mouse columns.');
    end

    % Infer sampling from the time column; size baseline / shuffle blocks to match
    timeCol = tbl{:, timeIdx};
    [timeMinutesAll, TsMinutes] = lb_infer_time_minutes(timeCol, varNames{timeIdx});
    opts.TsMinutes = TsMinutes;
    opts.blockLenSamples = max(4, round((opts.blockLenHours * 60) / TsMinutes));
    durationHours = (max(timeMinutesAll) - min(timeMinutesAll)) / 60;
    % Baseline window: roughly half the recording, clamped between 72 and 168 h
    opts.baselineWinSamples = max(20, round((min(max(72, 0.5*durationHours), 168) * 60) / TsMinutes));

    LightDur_h = median(tbl{:, lightIdx}, 'omitnan');
    if ~isfinite(LightDur_h), LightDur_h = NaN; end

    % Condition groups for scalogram averages (define once, reuse thereafter)
    if isempty(groupsTemplate)
        groups = lb_group_dialog(varNames, mouseIdx);
        groupsTemplate = groups;
    else
        groups = lb_resolve_groups(varNames, groupsTemplate);
        fprintf('  Reusing %d condition group(s).\n', numel(groups));
    end
    for g = 1:numel(groups)
        groups(g).colNames = varNames(groups(g).colIdx);
    end

    % Output layout: Reports/ + TimeSeries/Residual|Removed/ + Figures_Scalograms/
    [~, fileStem, ~] = fileparts(inputFile);
    runFolder = fullfile(outputRoot, fileStem);
    reportsFolder = fullfile(runFolder, 'Reports');
    tsResFolder = fullfile(runFolder, 'TimeSeries', 'Residual');
    tsRemFolder = fullfile(runFolder, 'TimeSeries', 'Removed');
    lb_ensure_dir(reportsFolder);
    lb_ensure_dir(tsResFolder);
    lb_ensure_dir(tsRemFolder);

    % Pre-allocate residual matrices for each period floor (Min360, Min60)
    nMice = numel(mouseIdx);
    keys = arrayfun(@(mp) sprintf('Min%d', mp), opts.minPeriodsMinutes, 'UniformOutput', false);
    res = struct();
    for ki = 1:numel(keys)
        res.(keys{ki}).sel = NaN(height(tbl), nMice);
        res.(keys{ki}).full = NaN(height(tbl), nMice);
        res.(keys{ki}).remSel = NaN(height(tbl), nMice);
        res.(keys{ki}).remFull = NaN(height(tbl), nMice);
    end

    anchorRows = cell(nMice, 8);
    recRows = cell(nMice, 7);
    timeOut = tbl{:, timeIdx};
    lightOut = tbl{:, lightIdx};

    % --- Per-mouse analysis ---
    for c = 1:nMice
        colName = varNames{mouseIdx(c)};
        fprintf('  Mouse %d/%d: %s\n', c, nMice, colName);
        x0 = tbl{:, mouseIdx(c)};
        if ~isnumeric(x0), x0 = str2double(string(x0)); end

        mouseOut = lb_hsub_process_mouse(x0, timeMinutesAll, opts, LightDur_h);
        mouseOut.colName = colName;
        rec = lb_hsub_recommendation(mouseOut, opts);

        % Rows collected for the summary workbook
        anchorRows(c, :) = {colName, mouseOut.anchorOK, mouseOut.P0_h, ...
            mouseOut.maxDeltaR2, mouseOut.pBlock, mouseOut.cyclesAtP0, mouseOut.note, rec.code};
        recRows(c, :) = {colName, rec.code, rec.mode, rec.minKey, mouseOut.anchorOK, rec.reason, mouseOut.P0_h};

        for ki = 1:numel(keys)
            k = keys{ki};
            if isfield(mouseOut.residuals, k)
                res.(k).sel(:, c) = mouseOut.residuals.(k).selective;
                res.(k).full(:, c) = mouseOut.residuals.(k).full;
                res.(k).remSel(:, c) = mouseOut.residuals.(k).removedSel;
                res.(k).remFull(:, c) = mouseOut.residuals.(k).removedFull;
            end
        end
    end

    % --- Write residual / removed workbooks (all arms) ---
    mouseNames = varNames(mouseIdx);
    for ki = 1:numel(keys)
        k = keys{ki};
        % Column names: Residual_Selective_<MouseID>, etc.
        % The across-photoperiod step looks up Residual_Selective_Min360_* specifically.
        lb_write_ts(tsResFolder, sprintf('Residual_Selective_%s_%s.xlsx', k, fileStem), ...
            timeOut, lightOut, mouseNames, res.(k).sel, 'Residual_Selective');
        lb_write_ts(tsResFolder, sprintf('Residual_FullLadder_%s_%s.xlsx', k, fileStem), ...
            timeOut, lightOut, mouseNames, res.(k).full, 'Residual_FullLadder');
        lb_write_ts(tsRemFolder, sprintf('Removed_Selective_%s_%s.xlsx', k, fileStem), ...
            timeOut, lightOut, mouseNames, res.(k).remSel, 'Removed_Selective');
        lb_write_ts(tsRemFolder, sprintf('Removed_FullLadder_%s_%s.xlsx', k, fileStem), ...
            timeOut, lightOut, mouseNames, res.(k).remFull, 'Removed_FullLadder');
    end

    % One "recommended" series per mouse: usually SEL_P360; PASS → raw pass-through
    recMatrix = NaN(height(tbl), nMice);
    for c = 1:nMice
        minKey = recRows{c, 4};
        mode = recRows{c, 3};
        anchorOkThis = recRows{c, 5};
        if islogical(anchorOkThis)
            okFlag = anchorOkThis;
        else
            okFlag = logical(anchorOkThis);
        end
        if ~okFlag
            recMatrix(:, c) = tbl{:, mouseIdx(c)};   % anchor failed → keep raw
        elseif strcmpi(mode, 'FullLadder')
            recMatrix(:, c) = res.(minKey).full(:, c);
        else
            recMatrix(:, c) = res.(minKey).sel(:, c);
        end
    end
    lb_write_ts(tsResFolder, sprintf('Residual_Recommended_%s.xlsx', fileStem), ...
        timeOut, lightOut, mouseNames, recMatrix, 'Residual_Recommended');

    % Anchor + recommendation sheets used later as validation metadata
    anchorTable = cell2table(anchorRows, 'VariableNames', ...
        {'Column', 'AnchorOK', 'Period_hours', 'MaxDeltaR2', 'pBlock', 'Cycles', 'Notes', 'Recommendation'});
    recTable = cell2table(recRows, 'VariableNames', ...
        {'Column', 'Code', 'Mode', 'MinKey', 'AnchorOK', 'Reason', 'Period_hours'});
    summaryXlsx = fullfile(reportsFolder, 'HarmonicRemoval_Summary.xlsx');
    if isfile(summaryXlsx), delete(summaryXlsx); end
    writetable(anchorTable, summaryXlsx, 'Sheet', 'Anchor_Report');
    writetable(recTable, summaryXlsx, 'Sheet', 'Recommendation');

    fprintf('  Wrote: %s\n', summaryXlsx);

    % After residuals exist on disk, draw SEL_P360 Removed|Residual scalograms.
    % Only AnchorOK mice enter the group average (failed anchors stay out).
    armKey = opts.scalogramArm;
    if isfield(res, armKey)
        anchorOK = false(nMice, 1);
        for c = 1:nMice
            if islogical(recRows{c, 5})
                anchorOK(c) = recRows{c, 5};
            else
                anchorOK(c) = logical(recRows{c, 5});
            end
        end
        try
            lb_hsub_plot_scalograms(tbl, timeIdx, lightIdx, varNames, groups, ...
                res.(armKey).sel, res.(armKey).remSel, mouseIdx, anchorOK, ...
                plotMode, saveIndividualScalograms, fileStem, runFolder, ...
                LightDur_h, opts);
        catch ME
            % Scalograms need the Wavelet Toolbox; residuals/reports still stand alone
            warning('HSub scalograms failed for %s: %s', fileStem, ME.message);
        end
    end
end

function out = lb_hsub_process_mouse(x0, timeMinutesAll, opts, LightDur_h)
% One mouse column: detrending, anchor scan, surrogate gate, selective/full residual.
% If the anchor fails, residuals stay as the original series (pass-through).

    out = struct('anchorOK', false, 'P0_h', NaN, 'maxDeltaR2', NaN, 'pBlock', NaN, ...
        'cyclesAtP0', NaN, 'peakZ', NaN, 'note', '', 'residuals', struct(), 'crossQC', struct());

    % --- Cleanup: optional short-gap fill, then slow baseline removal ---
    x0 = x0(:);
    missingOriginal = isnan(x0) | ~isfinite(x0);
    missingFrac = sum(missingOriginal) / max(numel(x0), 1);
    maxGapSamples = max(0, round(opts.maxGapMinutes / opts.TsMinutes));
    xFill = x0;
    if missingFrac <= opts.missingFracThreshold && maxGapSamples > 0
        xFill = fillmissing(x0, 'pchip', 'MaxGap', maxGapSamples);
    end

    % Moving median ≈ multi-day trend; subtract so sinusoid fits see the rhythm
    base = movmedian(xFill, opts.baselineWinSamples, 'omitnan', 'Endpoints', 'shrink');
    yDetrFull = xFill - base;
    validIdx = find(~isnan(yDetrFull) & isfinite(timeMinutesAll));

    % Default residual = raw (overwritten only if AnchorOK and subtraction runs)
    keys = arrayfun(@(mp) sprintf('Min%d', mp), opts.minPeriodsMinutes, 'UniformOutput', false);
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
    durCol_h = (max(t_min) - min(t_min)) / 60;
    y = yDetrFull(validIdx);
    y = y - median(y, 'omitnan');   % centre before ΔR² / sinusoid fits

    % --- Circadian anchor: P0 = period that maximises ΔR² in 22–28 h ---
    % ΔR² = R²(drift + sinusoid at P) − R²(drift only)
    periods_h = (opts.anchorBandHours(1):opts.anchorStepHours:opts.anchorBandHours(2));
    [P0_h, maxDeltaR2, peakZ] = lb_anchor_scan(t_hr, y, periods_h);
    cyclesAtP0 = durCol_h / max(P0_h, eps);

    blk = opts.blockLenSamples;
    if numel(y) < 4 * blk, blk = max(4, floor(numel(y) / 4)); end

    % --- Surrogate gate: shuffle blocks, re-scan max ΔR² under the null ---
    % Block shuffle keeps short-range correlation but breaks circadian phase.
    surrMax = NaN(opts.nSurrogates, 1);
    for s = 1:opts.nSurrogates
        yS = lb_block_shuffle(y, blk);
        surrMax(s) = lb_anchor_max_delta(t_hr, yS, periods_h);
    end
    % Monte Carlo p-value (add-one smoothing)
    pBlock = (1 + sum(surrMax >= maxDeltaR2)) / (1 + opts.nSurrogates);

    % Reject peaks that sit on the edge of the search band (ambiguous maximum)
    edgeOK = true;
    if opts.useEdgeGuard
        edgeOK = (P0_h > opts.anchorBandHours(1) + opts.edgeMarginHours) && ...
            (P0_h < opts.anchorBandHours(2) - opts.edgeMarginHours);
    end

    anchorOK = isfinite(P0_h) && isfinite(maxDeltaR2) && isfinite(pBlock) && ...
        (pBlock < opts.alphaAnchor) && (cyclesAtP0 >= opts.minCyclesForAnchor) && ...
        (maxDeltaR2 >= opts.minDeltaR2) && edgeOK;

    out.anchorOK = anchorOK;
    out.P0_h = P0_h;
    out.maxDeltaR2 = maxDeltaR2;
    out.peakZ = peakZ;
    out.pBlock = pBlock;
    out.cyclesAtP0 = cyclesAtP0;
    out.duration_h = durCol_h;
    out.LightDur_h = LightDur_h;

    if ~anchorOK
        % Residuals already set to raw above — nothing further to subtract
        out.note = sprintf('Anchor rejected (pBlock=%.3g, cycles=%.2f, dR2=%.4g)', ...
            pBlock, cyclesAtP0, maxDeltaR2);
        return;
    end

    % --- Sliding windows: track fundamental A, φ at P0 ---
    winHours = min(max(72, round(durCol_h / 2)), 168);
    [winStarts, winEnds] = lb_make_windows(min(t_min), max(t_min), winHours*60, opts.stepHours*60);
    nW = numel(winStarts);
    P0_min = P0_h * 60;

    A1_w = NaN(nW, 1); Phi1_w = NaN(nW, 1);
    for w = 1:nW
        wMask = (t_min >= winStarts(w)) & (t_min < winEnds(w));
        if sum(wMask) < opts.minSamplesPerWindow, continue; end
        [A1, phi1] = lb_fit_sinusoid(t_min(wMask), y(wMask), P0_min);
        A1_w(w) = A1; Phi1_w(w) = phi1;
    end

    % Long vs short recordings change how strictly we call a harmonic "likely"
    nFund = sum(isfinite(A1_w) & isfinite(Phi1_w));
    if nFund >= 4
        regime = 'long'; minWin = 4;
    elseif nFund >= 2
        regime = 'short'; minWin = 2;
    else
        regime = 'insufficient'; minWin = Inf;
    end

    % --- For each period floor (360 min, 60 min): screen harmonics, subtract ---
    for mp = opts.minPeriodsMinutes
        key = sprintf('Min%d', mp);
        % Highest harmonic order allowed by this period floor
        K = max(2, min(24, floor(P0_min / mp)));
        kList = 2:K;

        Ak_w = NaN(nW, numel(kList));
        Phik_w = NaN(nW, numel(kList));
        ratio_w = NaN(nW, numel(kList));

        % In each window, fit each harmonic k near P0/k (narrow period scan)
        for w = 1:nW
            wMask = (t_min >= winStarts(w)) & (t_min < winEnds(w));
            if sum(wMask) < opts.minSamplesPerWindow || ~isfinite(A1_w(w)), continue; end
            tmw = t_min(wMask); yw = y(wMask);
            for j = 1:numel(kList)
                k = kList(j);
                targetP = P0_min / k;
                scanP = linspace(targetP*(1-opts.peakSearchFrac), targetP*(1+opts.peakSearchFrac), 19);
                bestR2 = -Inf; bestA = NaN; bestPhi = NaN; bestP = NaN;
                for sp = 1:numel(scanP)
                    [Akk, phikk, R2kk] = lb_fit_sinusoid(tmw, yw, scanP(sp));
                    if isfinite(R2kk) && R2kk > bestR2
                        bestR2 = R2kk; bestA = Akk; bestPhi = phikk; bestP = scanP(sp);
                    end
                end
                if ~isfinite(bestR2), continue; end
                Ak_w(w,j) = bestA; Phik_w(w,j) = bestPhi; ratio_w(w,j) = bestP / targetP;
            end
        end

        % Decide which harmonics are "likely" (used by Selective arm)
        harmonicLikely = false(1, numel(kList));
        for j = 1:numel(kList)
            okW = isfinite(A1_w) & isfinite(Ak_w(:,j)) & isfinite(Phi1_w) & isfinite(Phik_w(:,j));
            if sum(okW) < minWin, continue; end
            % Period ratio ≈ 1; phase of k-th component locks to k × φ1
            ratioPass = mad(ratio_w(okW,j) - 1, 1) <= opts.ratioTol;
            del = lb_wrap_pi(Phik_w(okW,j) - kList(j) * Phi1_w(okW));
            Rplv = abs(mean(exp(1i * del)));   % phase-locking value
            Rpass = isfinite(Rplv) && (Rplv >= opts.Rcrit);
            if strcmp(regime, 'long')
                % Need ≥2 of: ratio, amplitude correlation, phase lock
                rho = corr(A1_w(okW), Ak_w(okW,j), 'Type', 'Spearman', 'Rows', 'complete');
                harmonicLikely(j) = (ratioPass + (isfinite(rho) && abs(rho) >= opts.rhoCrit) + Rpass) >= 2;
            else
                % Short records: require ratio + phase lock + relative amplitude
                A1_med = median(A1_w(okW), 'omitnan');
                Ak_med = median(Ak_w(okW,j), 'omitnan');
                ampPass = isfinite(A1_med) && A1_med > 0 && isfinite(Ak_med) && ...
                    (Ak_med >= opts.ampFracCritShort * A1_med);
                harmonicLikely(j) = ratioPass && Rpass && ampPass;
            end
        end

        kLikely = kList(harmonicLikely);
        if opts.subtractFundamentalInSelective
            kSubSel = unique([1, kLikely], 'stable');   % always include k=1
        else
            kSubSel = kLikely;
        end

        % Selective: least-squares fit of drift + chosen cos/sin columns; residual = y − harmonics
        if isempty(kSubSel)
            residualSel = y; removedSel = zeros(size(y)); varExplSel = 0;
        else
            [X0, Xh] = lb_build_design(t_hr, t_min, P0_min, kSubSel);
            beta = [X0, Xh] \ y;
            yhatDrift = X0 * beta(1:size(X0,2));
            yhatTotal = [X0, Xh] * beta;
            removedSel = yhatTotal - yhatDrift;   % harmonic part only (not the drift)
            residualSel = y - removedSel;
            varExplSel = 1 - sum((y-residualSel).^2,'omitnan') / (sum((y-mean(y,'omitnan')).^2,'omitnan')+eps);
        end

        % Full ladder: subtract k = 1…K regardless of the likelihood screen
        kSubFull = 1:K;
        [X0F, XhF] = lb_build_design(t_hr, t_min, P0_min, kSubFull);
        betaF = [X0F, XhF] \ y;
        yhatDriftF = X0F * betaF(1:size(X0F,2));
        yhatTotalF = [X0F, XhF] * betaF;
        removedFull = yhatTotalF - yhatDriftF;
        residualFull = y - removedFull;
        varExplFull = 1 - sum((y-residualFull).^2,'omitnan') / (sum((y-mean(y,'omitnan')).^2,'omitnan')+eps);

        % Map residuals back onto the full-length (including missing) series
        out.residuals.(key).selective = lb_map_full(x0, validIdx, residualSel);
        out.residuals.(key).full = lb_map_full(x0, validIdx, residualFull);
        out.residuals.(key).removedSel = lb_map_full(NaN(size(x0)), validIdx, removedSel);
        out.residuals.(key).removedFull = lb_map_full(NaN(size(x0)), validIdx, removedFull);

        % How different are selective vs full? Used later for the recommendation
        rmsDiff = sqrt(mean((residualSel - residualFull).^2, 'omitnan'));
        nRMS = rmsDiff / (sqrt(mean(y.^2, 'omitnan')) + eps);
        out.crossQC.(key).nRMS = nRMS;
        out.crossQC.(key).dVar = varExplFull - varExplSel;
    end

    out.note = sprintf('OK (%s regime)', regime);
end

function rec = lb_hsub_recommendation(mouseOut, opts)
% Choose which residual to recommend for this mouse.
% Default when AnchorOK: Selective + Min360 → code 'SEL_P360'.
% Tip toward Min60 or FullLadder only if selective/full diverge strongly.
    rec = struct('code', 'PASS', 'mode', 'Selective', 'minKey', 'Min360', 'reason', 'Anchor rejected');
    if ~mouseOut.anchorOK, return; end   % no reliable circadian → keep raw (PASS)

    % Prefer the 6 h period floor unless Min60 selective/full differ a lot
    chosen = 'Min360';
    if isfield(mouseOut, 'crossQC') && isfield(mouseOut.crossQC, 'Min60')
        if isfinite(mouseOut.crossQC.Min60.nRMS) && mouseOut.crossQC.Min60.nRMS >= opts.nRMS_high
            chosen = 'Min60';
        end
    end

    % Prefer Selective; switch to FullLadder if QC says the two arms disagree
    mode = 'Selective';
    if isfield(mouseOut.crossQC, chosen)
        qc = mouseOut.crossQC.(chosen);
        if qc.nRMS >= opts.nRMS_high || qc.dVar >= opts.dVar_high
            mode = 'FullLadder';
        end
    end

    rec.minKey = chosen;
    rec.mode = mode;
    % Build short code string, e.g. SEL_P360 or FULL_P60
    if strcmp(mode, 'FullLadder')
        rec.code = sprintf('FULL_%s', strrep(chosen, 'Min', 'P'));
    else
        rec.code = sprintf('SEL_%s', strrep(chosen, 'Min', 'P'));
    end
    rec.reason = sprintf('%s chosen', rec.code);
end

%% --------------------- maths helpers ------------------------------------

function [P0_h, maxDeltaR2, peakZ] = lb_anchor_scan(t_hr, y, periods_h)
% Try every candidate period in periods_h.
% Keep the period with the largest ΔR² (= best improvement from adding a sinusoid).
    delta = NaN(size(periods_h));
    for i = 1:numel(periods_h)
        delta(i) = lb_delta_r2_one(t_hr, y, periods_h(i));
    end
    [maxDeltaR2, ix] = max(delta);   % ix = index of the winning period
    P0_h = periods_h(ix);

    % How unusual is that peak within the scan? (descriptive QC only)
    mu = mean(delta, 'omitnan'); sd = std(delta, 0, 'omitnan');
    if isfinite(sd) && sd > 0
        peakZ = (maxDeltaR2 - mu) / sd;
    else
        peakZ = NaN;
    end
end

function d = lb_anchor_max_delta(t_hr, y, periods_h)
% Same scan as lb_anchor_scan, but only returns the max ΔR².
% Used inside the surrogate loop (faster / simpler — no need for P0).
    d = -Inf;
    for i = 1:numel(periods_h)
        d = max(d, lb_delta_r2_one(t_hr, y, periods_h(i)));
    end
end

function dR2 = lb_delta_r2_one(t_hr, y, P_h)
% Compare two least-squares fits at one candidate period P_h:
%   model 0 = slow drift only
%   model 1 = drift + cos/sin at period P_h
% ΔR² = R²_1 − R²_0  (how much better the sinusoid makes the fit).
    t_min = t_hr(:) * 60;
    y = y(:);
    tz = lb_zscore(t_hr(:));                          % standardised time (for drift)
    X0 = [ones(size(tz)), tz, tz.^2];                 % columns: intercept, time, time²
    w = 2 * pi / (P_h * 60);                          % angular frequency for P_h
    X1 = [X0, cos(w*t_min), sin(w*t_min)];            % drift + sinusoid

    % "\" = MATLAB least-squares solve for coefficients
    b0 = X0 \ y; b1 = X1 \ y;
    r0 = y - X0*b0; r1 = y - X1*b1;                   % residuals of each model
    ssTot = sum((y - mean(y,'omitnan')).^2, 'omitnan') + eps;
    R2_0 = 1 - sum(r0.^2, 'omitnan') / ssTot;
    R2_1 = 1 - sum(r1.^2, 'omitnan') / ssTot;
    dR2 = R2_1 - R2_0;
end

function [A, phi, R2] = lb_fit_sinusoid(t_min, y, period_min)
% Fit: y ≈ A*cos(2πt/P + φ) plus a simple linear drift.
% Returns amplitude A, phase phi, and R² of the fit.
    if nargout < 3, R2 = NaN; end
    t = t_min(:); y = y(:);
    if numel(t) < 10 || all(~isfinite(y)), A = NaN; phi = NaN; R2 = NaN; return; end
    w = 2 * pi / period_min;
    tZ = lb_zscore(t);
    X = [cos(w*t), sin(w*t), ones(size(t)), tZ];   % cos, sin, intercept, drift
    b = X \ y; yhat = X * b;
    % Convert cos/sin coefficients into amplitude and phase
    A = hypot(b(1), b(2));
    phi = atan2(-b(2), b(1));
    ssRes = sum((y - yhat).^2, 'omitnan');
    ssTot = sum((y - mean(y,'omitnan')).^2, 'omitnan');
    R2 = 1 - ssRes / (ssTot + eps);
end

function [X0, Xh] = lb_build_design(t_hr, t_min, P0_min, kList)
% Build the design matrix used when subtracting harmonics.
%   X0 = quadratic drift columns
%   Xh = cos/sin pair for each harmonic order in kList (k=1 is the fundamental)
    tz = lb_zscore(t_hr(:));
    X0 = [ones(size(tz)), tz, tz.^2];
    t = t_min(:);
    Xh = zeros(numel(t), 2*numel(kList));   % 2 columns per harmonic
    for ii = 1:numel(kList)
        w = 2 * pi * kList(ii) / P0_min;    % frequency of harmonic k
        col = 2*(ii-1) + 1;
        Xh(:, col) = cos(w*t);
        Xh(:, col+1) = sin(w*t);
    end
end

function yS = lb_block_shuffle(y, blockLen)
% Surrogate: cut the series into chunks of length blockLen, then shuffle
% the order of those chunks. Local wiggles stay; long circadian phase is broken.
    y = y(:); n = numel(y);
    if n < 2*blockLen
        yS = y(randperm(n));   % too short for blocks → shuffle every sample
        return;
    end
    nB = floor(n / blockLen); remN = n - nB*blockLen;
    blocks = cell(nB + (remN > 0), 1);
    for b = 1:nB
        blocks{b} = y((b-1)*blockLen+1 : b*blockLen);
    end
    if remN > 0, blocks{end} = y(nB*blockLen+1:end); end   % leftover tail
    yS = vertcat(blocks{randperm(numel(blocks))});         % reassemble in random order
end

function [starts, ends] = lb_make_windows(tMin, tMax, winLenMin, stepMin)
% Sliding window start/end times (in minutes) along the recording.
    starts = (tMin:stepMin:(tMax - winLenMin));
    ends = starts + winLenMin;
    if isempty(starts), starts = tMin; ends = tMax; end
end

function x = lb_wrap_pi(x)
% Force angles into (−π, π] so phase differences are comparable.
    x = mod(x + pi, 2*pi) - pi;
end

function z = lb_zscore(x)
% Standardise: (x − mean) / sd. If sd is 0, just subtract the mean.
    x = x(:); sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        z = x - mean(x, 'omitnan');
    else
        z = (x - mean(x, 'omitnan')) / sd;
    end
end

function full = lb_map_full(template, validIdx, shortVec)
% Analysis often drops NaN rows. This puts the short result vector back
% into a full-length series (same length as the original mouse column).
    full = template;
    if all(isnan(template)), full = NaN(size(template)); end
    full(validIdx) = shortVec;
end

%% --------------------- I/O helpers --------------------------------------

function [tbl, varNames] = lb_read_excel(inputFile)
% Read Excel/CSV into a MATLAB table; keep original header names.
    tbl = readtable(inputFile, 'VariableNamingRule', 'preserve');
    varNames = tbl.Properties.VariableNames;
end

function [timeIdx, lightIdx, mouseIdx] = lb_map_columns(varNames)
% Work out which columns are Time, Light duration, and mouse activity.
% mouseIdx may be refined by a list dialog if several candidates remain.
    timeIdx = []; lightIdx = [];
    for i = 1:numel(varNames)
        vn = lower(strtrim(varNames{i}));
        if contains(vn, 'time') && isempty(timeIdx), timeIdx = i; end
        if contains(vn, 'light') && contains(vn, 'duration'), lightIdx = i; end
    end
    % Everything else is a candidate mouse column
    mouseIdx = setdiff(1:numel(varNames), [timeIdx, lightIdx]);
    keep = true(size(mouseIdx));
    for k = 1:numel(mouseIdx)
        vn = lower(varNames{mouseIdx(k)});
        % Drop columns that are clearly not activity
        if contains(vn, 'zt') || strcmp(vn, 'day') || contains(vn, 'unnamed')
            keep(k) = false;
        end
    end
    mouseIdx = mouseIdx(keep);
    if numel(mouseIdx) > 1
        [sel, ok] = listdlg('PromptString', 'Select mouse activity columns:', ...
            'SelectionMode', 'multiple', 'ListString', varNames(mouseIdx), ...
            'InitialValue', 1:numel(mouseIdx));
        if ok && ~isempty(sel), mouseIdx = mouseIdx(sel); end
    end
end

function [timeMinutes, TsMinutes] = lb_infer_time_minutes(timeCol, timeName)
% Convert the time column to minutes and estimate the sampling step.
% If the header says "hr" (or values look like hours), multiply by 60.
    if ~isnumeric(timeCol), timeCol = str2double(string(timeCol)); end
    timeCol = timeCol(:);
    if contains(lower(timeName), 'hr') || (max(timeCol, [], 'omitnan') < 500)
        timeMinutes = timeCol * 60;
    else
        timeMinutes = timeCol;
    end
    d = diff(timeMinutes);
    d = d(isfinite(d) & d > 0);
    if isempty(d)
        TsMinutes = 10;          % fallback: assume 10 min sampling
    else
        TsMinutes = median(d);   % typical step between samples
    end
end

function lb_write_ts(folder, fileName, timeOut, lightOut, mouseNames, dataMatrix, prefix)
% Write one timeseries Excel file:
%   Time | prefix_Mouse1 | prefix_Mouse2 | … | Light_duration_h
    if nargin < 7 || isempty(prefix), prefix = 'Value'; end
    headers = [{'Time'}, strcat(prefix, '_', mouseNames), {'Light_duration_h'}];
    T = array2table([timeOut, dataMatrix, lightOut], ...
        'VariableNames', matlab.lang.makeValidName(headers));
    writetable(T, fullfile(folder, fileName));
end

function lb_ensure_dir(p)
% Create folder p if it does not already exist.
    if ~exist(p, 'dir'), mkdir(p); end
end

%% --------------------- Groups + HSub scalograms -------------------------
% Condition groups (e.g. F / M) drive the group-average Removed|Residual
% scalograms saved under Figures_Scalograms/. The same group names are
% reused across a multi-file batch so photoperiods stay aligned.

function groups = lb_group_dialog(varNames, mouseIdx)
% Interactive condition groups (same flow as the wavelet code):
%   1. Type how many groups
%   2. For each group: type the name, then tick which mouse columns belong to it
%
% Returns a struct array with fields: name, colIdx, colNames.
% If the user cancels the count dialog, everything is put in one group "All".
    nDef = '2';
    ansN = inputdlg({'Number of condition groups (e.g. 2 for F and M):'}, ...
        'Groups', 1, {nDef});
    if isempty(ansN)
        groups = struct('name', {'All'}, 'colIdx', {mouseIdx}, 'colNames', {varNames(mouseIdx)});
        return;
    end
    nG = max(1, round(str2double(ansN{1})));
    if ~isfinite(nG), nG = 1; end

    groups = struct('name', {}, 'colIdx', {}, 'colNames', {});
    remaining = mouseIdx;   % columns not yet assigned to a group
    for g = 1:nG
        nameAns = inputdlg({sprintf('Name for group %d/%d:', g, nG)}, ...
            'Group name', 1, {sprintf('Group%d', g)});
        if isempty(nameAns)
            grpName = sprintf('Group%d', g);
        else
            grpName = strtrim(nameAns{1});
            if isempty(grpName), grpName = sprintf('Group%d', g); end
        end
        if isempty(remaining)
            warning('No columns left for group "%s".', grpName);
            break;
        end
        [sel, ok] = listdlg('PromptString', sprintf('Select columns for "%s":', grpName), ...
            'SelectionMode', 'multiple', 'ListString', varNames(remaining), ...
            'InitialValue', 1:numel(remaining));
        if ~ok || isempty(sel)
            warning('No columns selected for group "%s"; skipped.', grpName);
            continue;
        end
        colIdx = remaining(sel);
        groups(end+1).name = grpName; %#ok<AGROW>
        groups(end).colIdx = colIdx;
        groups(end).colNames = varNames(colIdx);
        remaining = setdiff(remaining, colIdx, 'stable');
    end
    if isempty(groups)
        groups = struct('name', {'All'}, 'colIdx', {mouseIdx}, 'colNames', {varNames(mouseIdx)});
    end
end

function groups = lb_resolve_groups(varNames, template)
% Re-map template groups onto this file's column names.
% Matching is by mouse column *name* (not column index), so files with the
% same animals in a different order still line up.
    groups = struct('name', {}, 'colIdx', {}, 'colNames', {});
    for g = 1:numel(template)
        names = template(g).colNames;
        if isstring(names), names = cellstr(names); end
        colIdx = [];
        for k = 1:numel(names)
            ix = find(strcmp(varNames, names{k}), 1);
            if ~isempty(ix), colIdx(end+1) = ix; end %#ok<AGROW>
        end
        groups(g).name = template(g).name;
        groups(g).colIdx = colIdx;
        groups(g).colNames = varNames(colIdx);
    end
end

function figurePaths = lb_hsub_plot_scalograms(tbl, timeIdx, lightIdx, varNames, groups, ...
        residualMat, removedMat, mouseIdx, anchorOK, plotMode, saveIndividual, ...
        fileStem, runFolder, LightDur_h, opts)
% Two-panel Removed|Residual scalograms for the selective Min360 arm.
%
% Left panel  = fitted circadian + harmonics that were subtracted ("Removed")
% Right panel = what remains after subtraction ("Residual")
%
% Group averages use AnchorOK mice only. Optional individuals are written
% under Group_*/Individuals/ at a slightly higher DPI for QC.
    figurePaths = {};
    theme = lb_hsub_theme(plotMode, opts);

    % Time axis in days; white dotted lines mark light-duration changes if any
    time_hours = tbl{:, timeIdx};
    if ~isnumeric(time_hours), time_hours = str2double(string(time_hours)); end
    time_min = time_hours * 60;
    time_day = time_min / (60 * 24);
    TsMinutes = median(diff(time_min), 'omitnan');
    if ~isfinite(TsMinutes) || TsMinutes <= 0, TsMinutes = opts.TsMinutes; end
    lightVec = tbl{:, lightIdx};
    if ~isnumeric(lightVec), lightVec = str2double(string(lightVec)); end
    condChangeIdx = find(diff(lightVec) ~= 0);

    scaloRoot = fullfile(runFolder, 'Figures_Scalograms');
    lb_ensure_dir(scaloRoot);
    armLabel = opts.scalogramLabel;
    photoTag = sprintf('L%g', LightDur_h);
    if ~isfinite(LightDur_h), photoTag = fileStem; end

    FB = lb_hsub_make_filterbank(height(tbl), TsMinutes, opts);
    fprintf('  HSub scalograms (%s):\n', armLabel);

    for g = 1:numel(groups)
        if isempty(groups(g).colIdx), continue; end
        grpName = groups(g).name;

        % Mask: mouse is in this group AND passed the circadian anchor gate
        colMask = false(1, numel(mouseIdx));
        for c = 1:numel(mouseIdx)
            if ~anchorOK(c), continue; end
            if any(mouseIdx(c) == groups(g).colIdx)
                colMask(c) = true;
            end
        end
        nOk = sum(colMask);
        if nOk == 0
            warning('Group "%s": no AnchorOK mice; skipped HSub scalogram.', grpName);
            continue;
        end

        grpFolder = fullfile(scaloRoot, ['Group_' lb_sanitise(grpName)]);
        lb_ensure_dir(grpFolder);

        % Mean across AnchorOK mice, then one CWT each for Removed and Residual
        avgRemoved = mean(removedMat(:, colMask), 2);
        avgResidual = mean(residualMat(:, colMask), 2);
        [wtRem, periods_hours] = lb_hsub_cwt(avgRemoved, FB);
        [wtRes, ~] = lb_hsub_cwt(avgResidual, FB);

        outFile = fullfile(grpFolder, sprintf('HSub_Removed_Residual_%s_%s_%s.%s', ...
            armLabel, lb_sanitise(grpName), photoTag, theme.format));
        titleStr = sprintf('HSub %s | Average | %s | n=%d AnchorOK | L=%g h', ...
            armLabel, grpName, nOk, LightDur_h);
        lb_hsub_two_panel_scalogram(wtRem, wtRes, periods_hours, time_day, ...
            condChangeIdx, titleStr, outFile, theme, opts);
        figurePaths{end+1} = outFile; %#ok<AGROW>
        fprintf('    "%s" (n=%d) → %s\n', grpName, nOk, outFile);

        if saveIndividual
            indivFolder = fullfile(grpFolder, 'Individuals');
            lb_ensure_dir(indivFolder);
            indivTheme = theme;
            indivTheme.dpi = opts.indivScalogram.dpi;
            indivTheme.format = opts.indivScalogram.format;
            for c = find(colMask)
                mouseName = varNames{mouseIdx(c)};
                [wtRemI, ~] = lb_hsub_cwt(removedMat(:, c), FB);
                [wtResI, ~] = lb_hsub_cwt(residualMat(:, c), FB);
                outInd = fullfile(indivFolder, sprintf('HSub_Removed_Residual_%s_%s_%s.%s', ...
                    armLabel, lb_sanitise(grpName), lb_sanitise(mouseName), indivTheme.format));
                titleInd = sprintf('HSub %s | %s | %s | L=%g h', ...
                    armLabel, grpName, mouseName, LightDur_h);
                lb_hsub_two_panel_scalogram(wtRemI, wtResI, periods_hours, time_day, ...
                    condChangeIdx, titleInd, outInd, indivTheme, opts);
                figurePaths{end+1} = outInd; %#ok<AGROW>
            end
        end
    end
end

function theme = lb_hsub_theme(plotMode, opts)
% Pick export format/DPI: development = quick PNG; publication = high-res JPEG.
    if strcmpi(plotMode, 'publication')
        theme = opts.plot.publication;
    else
        theme = opts.plot.development;
    end
    theme.colormap = opts.colormap;
end

function FB = lb_hsub_make_filterbank(nSamples, TsMinutes, opts)
% Morlet filter bank for this recording length / sampling rate.
    FB = cwtfilterbank('SignalLength', nSamples, ...
        'SamplingPeriod', minutes(TsMinutes), ...
        'Wavelet', opts.waveletName, ...
        'PeriodLimits', minutes(opts.periodLimitsMinutes), ...
        'VoicesPerOctave', 32);
end

function [wt, periods_hours] = lb_hsub_cwt(signal, FB)
% Continuous wavelet transform; non-finite samples zeroed so cwt does not fail.
    signal = signal(:);
    signal(~isfinite(signal)) = 0;
    [wt, periods] = cwt(signal, 'FilterBank', FB);
    periods_hours = hours(periods);
end

function lb_hsub_two_panel_scalogram(wtRem, wtRes, periods_hours, time_day, ...
        condChangeIdx, titleStr, outPath, theme, opts)
% Side-by-side jet images: Removed (left) | Residual (right).
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1600 640]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    yLo = min(periods_hours); yHi = max(periods_hours);
    lb_hsub_draw_panel(nexttile(tl, 1), wtRem, periods_hours, time_day, ...
        condChangeIdx, yLo, yHi, 'Removed', theme, opts);
    lb_hsub_draw_panel(nexttile(tl, 2), wtRes, periods_hours, time_day, ...
        condChangeIdx, yLo, yHi, 'Residual', theme, opts);
    sgtitle(fig, titleStr, 'Interpreter', 'none', 'FontName', 'Arial');
    lb_hsub_export_fig(fig, outPath, theme);
    close(fig);
end

function lb_hsub_draw_panel(ax, wt, periods_hours, time_day, condChangeIdx, ...
        yLo, yHi, panelLabel, theme, opts)
% One jet |wt| panel; white dotted lines mark light-duration changes if present.
    axes(ax); %#ok<LAXES>
    pcolor(time_day, periods_hours, abs(wt));
    shading interp;
    colormap(ax, theme.colormap);
    colorbar;
    set(ax, 'YTick', opts.scalogramYTicksHours, 'FontName', 'Arial', 'FontSize', 11);
    hold(ax, 'on');
    for k = 1:numel(condChangeIdx)
        row = condChangeIdx(k) + 1;
        if row >= 1 && row <= numel(time_day)
            plot(ax, [time_day(row) time_day(row)], [yLo yHi], 'w:', 'LineWidth', 1.5);
        end
    end
    hold(ax, 'off');
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Period (hr)');
    title(ax, panelLabel);
end

function lb_hsub_export_fig(fig, outPath, theme)
% Save figure as PNG or JPEG at the chosen DPI.
    [folder, ~, ~] = fileparts(outPath);
    if ~isempty(folder), lb_ensure_dir(folder); end
    fmt = lower(theme.format);
    if strcmp(fmt, 'jpg'), fmt = 'jpeg'; end
    print(fig, outPath, ['-d' fmt], sprintf('-r%d', theme.dpi));
end

function s = lb_sanitise(strIn)
% Make a string safe for use in filenames (replace odd characters with _).
    s = char(regexprep(string(strIn), '[^\w\-]', '_'));
end
