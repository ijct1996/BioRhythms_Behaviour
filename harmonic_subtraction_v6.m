% -------------------------------------------------------------------------
% Script: Harmonic Subtration_v6
% -------------------------------------------------------------------------

clearvars; close all; clc;

%% ----------------------------- USER INPUTS ------------------------------

[fileName, filePath] = uigetfile({'*.xlsx','Excel files (*.xlsx)'}, 'Select the input .xlsx file');
if isequal(fileName, 0)
    fprintf('No input file selected. Exiting.\n');
    return;
end
inputFile = fullfile(filePath, fileName);
fprintf('Input file selected: %s\n', inputFile);

outputFolder = uigetdir(pwd, 'Select (or create) the output folder');
if isequal(outputFolder, 0)
    fprintf('No output folder selected. Exiting.\n');
    return;
end
fprintf('Output folder selected: %s\n', outputFolder);

%% ------------------------------ SUBFOLDERS ------------------------------

reportsFolder    = fullfile(outputFolder, 'Reports');

timeSeriesFolder = fullfile(outputFolder, 'TimeSeries');
tsResidualFolder = fullfile(timeSeriesFolder, 'Residual');
tsRemovedFolder  = fullfile(timeSeriesFolder, 'Removed');

figAnchorFolder  = fullfile(outputFolder, 'Figures_Anchor');
figScaloFolder   = fullfile(outputFolder, 'Figures_Scalograms');
figScaloResFolder = fullfile(figScaloFolder, 'Residual');
figScaloRemFolder = fullfile(figScaloFolder, 'Removed');

if ~exist(reportsFolder, 'dir');    mkdir(reportsFolder);    end
if ~exist(timeSeriesFolder, 'dir'); mkdir(timeSeriesFolder); end
if ~exist(tsResidualFolder, 'dir'); mkdir(tsResidualFolder); end
if ~exist(tsRemovedFolder,  'dir'); mkdir(tsRemovedFolder);  end
if ~exist(figAnchorFolder, 'dir');  mkdir(figAnchorFolder);  end

% scalograms folders are created only if enabled later

summaryXLSX = fullfile(reportsFolder, 'HarmonicRemoval_Reports_Summary.xlsx');
detailXLSX  = fullfile(reportsFolder, 'HarmonicRemoval_Reports_Detail.xlsx');

% Reset workbooks each run (clean, predictable)
if exist(summaryXLSX, 'file'); delete(summaryXLSX); end
if exist(detailXLSX,  'file'); delete(detailXLSX);  end

%% ------------------------- READ DATA (PRESERVE) --------------------------

try
    opts = detectImportOptions(inputFile, 'FileType', 'spreadsheet');
    opts.VariableNamingRule = 'preserve';
    opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA', '', 'N/A'});
    opts.ImportErrorRule = 'omitrow';
    opts.MissingRule = 'omitrow';
    dataTable = readtable(inputFile, opts);
catch ME
    fprintf('Error reading file: %s\n', ME.message);
    return;
end

if isempty(dataTable) || width(dataTable) < 3
    fprintf('Insufficient columns (need Time, >=1 Data, Light duration condition).\n');
    return;
end

varNames = dataTable.Properties.VariableNames;

%% -------------------------- COLUMN SELECTION UI --------------------------

[timeIdx, okTime] = listdlg('PromptString','Select the Time column', ...
    'SelectionMode','single','ListString',varNames,'InitialValue',1);
if ~okTime, fprintf('No Time column selected. Exiting.\n'); return; end

[lightIdx, okLight] = listdlg('PromptString','Select the Light duration condition column', ...
    'SelectionMode','single','ListString',varNames,'InitialValue',width(dataTable));
if ~okLight, fprintf('No Light duration column selected. Exiting.\n'); return; end

defaultData = setdiff(1:width(dataTable), [timeIdx, lightIdx], 'stable');
[dataIdx, okData] = listdlg('PromptString','Select the Data columns ', ...
    'SelectionMode','multiple','ListString',varNames,'InitialValue',defaultData);
if ~okData || isempty(dataIdx), fprintf('No Data columns selected. Exiting.\n'); return; end

dataIdx = setdiff(dataIdx, [timeIdx, lightIdx], 'stable');
if isempty(dataIdx), fprintf('After excluding Time and Light condition, no Data columns remain. Exiting.\n'); return; end

% -------------------------- OPTIONAL EXCLUDE STEP --------------------------

selectedDataNames = varNames(dataIdx);
excludeList = [{'<<None>>'}; selectedDataNames(:)];

[excludeSel, okExclude] = listdlg('PromptString', ...
    sprintf('Select Data columns to exclude.\nKeep "<<None>>" selected to exclude none.'), ...
    'SelectionMode','multiple', ...
    'ListString', excludeList, ...
    'InitialValue', 1);

excludedNames = {};
if okExclude
    excludeSel = setdiff(excludeSel, 1); % remove "<<None>>" if selected
    if ~isempty(excludeSel)
        excludePos = excludeSel - 1; % offset for "<<None>>"
        excludePos = excludePos(excludePos >= 1 & excludePos <= numel(dataIdx));
        if ~isempty(excludePos)
            excludedNames = varNames(dataIdx(excludePos));
            dataIdx(excludePos) = [];
        end
    end
end

if isempty(dataIdx)
    fprintf('All Data columns were excluded. Exiting.\n');
    return;
end

timeName  = varNames{timeIdx};
lightName = varNames{lightIdx};

fprintf('\nSelected columns:\n');
fprintf('  Time : %s\n', timeName);
fprintf('  Light duration condition (h): %s\n', lightName);
fprintf('  Data : %d columns\n', numel(dataIdx));
if ~isempty(excludedNames)
    fprintf('  Excluded : %d columns\n', numel(excludedNames));
end

%% ------------------------- FIXED INTERNAL SPECS --------------------------

% Always run both minimum periods (minutes)
minPeriodsToRun = [360, 60];   % Min360 then Min60

% Missingness rule
missingFracThreshold = 0.05;   % 5%
MaxGapMinutes        = 30;     % interpolate short gaps only (if <= 5%)

% Anchor scan band and resolution
anchorBand_h  = [22, 28];
anchorStep_h  = 0.01;          % 0.6 minutes

% Block-shuffle surrogate settings
blockLenHours = 8;
Nsurr         = 300;
alphaAnchor   = 0.05;

% Anchor viability (hard guards only)
minCyclesForAnchor = 3.5;
minDeltaR2         = 0.005;
useEdgeGuard       = true;
edgeMargin_h       = 0.10;

% QC flags (NOT hard rejection)
PeakZ_warn = 1.0;
SNR_warn   = 1.0;

% Harmonic coupling settings (Selective)
peakSearchFrac = 0.06;     % +/- 6% around target harmonic period for scan
ratioTol       = 0.08;     % MAD tolerance for ratio consistency
Rcrit          = 0.55;     % phase-locking criterion
rhoCrit        = 0.45;     % amplitude coupling (long regime)
ampFracCritShort = 0.04;   % short regime: Ak_med >= 4% of A1_med

% Windowing (dynamic)
stepHours = 24;
minSamplesForAnalysis = 250;
minSamplesPerWindow   = 80;

% Selective subtraction behaviour
SUBTRACT_FUNDAMENTAL_IN_SELECTIVE = true;

% Outputs: figures
SAVE_ANCHOR_FIGS = true;

% Outputs: scalograms
SAVE_SCALOGRAMS = true;
SAVE_SCALOGRAMS_SELECTIVE = true;  % ON: mirror FullLadder scalograms for Selective too

% Scalogram settings (match your wavelet script style)
SCALO_MIN_PERIOD_MIN = 60;    % minutes
SCALO_MAX_PERIOD_MIN = 1590;  % minutes (26.5 h)
SCALO_YTICKS_H = 0:4:26;

% QC thresholds (used for labels only; no popups)
nRMS_low  = 0.10;
nRMS_high = 0.25;
dVar_low  = 0.05;
dVar_high = 0.15;

%% --------------------- INFER TIME BASE AND SAMPLING ----------------------

timeCol = dataTable{:, timeIdx};
try
    [timeMinutesAll, TsMinutes] = infer_time_minutes(timeCol, timeName);
catch ME
    fprintf('Error inferring time base: %s\n', ME.message);
    return;
end

if ~isfinite(TsMinutes) || TsMinutes <= 0
    fprintf('Invalid sampling interval. Exiting.\n');
    return;
end

N = height(dataTable);
maxGapSamples = max(0, round(MaxGapMinutes / TsMinutes));

durationMinAll = max(timeMinutesAll, [], 'omitnan') - min(timeMinutesAll, [], 'omitnan');
durationHoursAll = durationMinAll / 60;

baselineWinHours = min(max(72, 0.5*durationHoursAll), 168);
baselineMinPoints = 20;
baselineWinSamples = max(baselineMinPoints, round((baselineWinHours*60) / TsMinutes));

blockLenSamples = max(4, round((blockLenHours*60) / TsMinutes));

fprintf('\nInferred sampling interval: %.6g minutes per sample\n', TsMinutes);
fprintf('Recording duration (approx): %.2f hours\n', durationHoursAll);
fprintf('Baseline removal: moving median window = %.1f h (%d samples)\n', baselineWinHours, baselineWinSamples);
fprintf('Anchor band fixed: %.1f-%.1f h (step %.3f h)\n', anchorBand_h(1), anchorBand_h(2), anchorStep_h);
fprintf('Block-shuffle surrogates: %d surrogates | block = %g h (%d samples)\n', Nsurr, blockLenHours, blockLenSamples);
fprintf('FullLadder runs: MinPeriods = [%s] minutes\n', num2str(minPeriodsToRun));

%% ------------------------- LIGHT CONDITION CHANGES -----------------------

% For scalogram vertical markers: find changes in Light duration condition (h)
lightVec = dataTable{:, lightIdx};
condChangeIdx = [];
try
    if isnumeric(lightVec) || islogical(lightVec)
        condChangeIdx = find(diff(lightVec) ~= 0);
    else
        s = string(lightVec);
        condChangeIdx = find(s(2:end) ~= s(1:end-1));
    end
catch
    condChangeIdx = [];
end

% time vectors for plotting scalograms
time_min = timeMinutesAll(:);
time_day = time_min / (60*24);

%% ---------------------------- REPORT TABLES ------------------------------

nCols = numel(dataIdx);

gapReportHeader = {'LightDur_h','Column','MissingFraction','Interpolated','FilledSamples','FilledGapsCount','LongGapsCount','UnfilledSamples','LongestGapSamples'};
gapReport = cell(nCols, numel(gapReportHeader));
gap_i = 0;

anchorHeader = {'LightDur_h','Column','AnchorOK','Period_hours','MaxDeltaR2','PeakZ','AnchorAmp','AmpSNR','Duration_h','CyclesAtPeriod','pBlock','Notes'};
anchorReport = cell(nCols, numel(anchorHeader));
anchor_i = 0;

excludedHeader = {'Column','Reason'};
excludedReport = cell(max(1, numel(excludedNames)), numel(excludedHeader));
for iE = 1:numel(excludedNames)
    excludedReport(iE,:) = {excludedNames{iE}, 'User excluded'};
end
if isempty(excludedNames)
    excludedReport = {'(none)','(none)'};
end

% Per-minPeriod reports
selHeader  = {'LightDur_h','MinPeriod_mins','Column','AnchorOK','Period_hours','RemovedFundamental','K_removed','Harmonics_k','VarExplained','NumWindowsUsed','Notes'};
fullHeader = {'LightDur_h','MinPeriod_mins','Column','AnchorOK','Period_hours','K_removed','VarExplained','Notes'};

selReport = struct();
fullReport = struct();
selWinReport = struct();

selWinHeader = {'LightDur_h','MinPeriod_mins','Column','WindowIndex','WinStart_h','WinEnd_h','Period_hours','k','Amp','Phase_rad','FitR2','PkBest_h','Ratio','Note'};

% Cross-QC
qcSelFullHeader = {'LightDur_h','MinPeriod_mins','Column','AnchorOK','Period_hours', ...
    'VarExpl_Sel','VarExpl_Full','DeltaVarExpl_FullMinusSel','RMSdiff','nRMS','Risk','Notes'};
qcSelFull = struct();

crossMinHeader = {'LightDur_h','Column','AnchorOK','Period_hours', ...
    'K_Min360','K_Min60','VarExplFull_Min360','VarExplFull_Min60','DeltaVarExpl_60minus360','RMS_60vs360','nRMS_60vs360','Recommendation','Notes'};
crossMinQC = cell(nCols, numel(crossMinHeader));
cross_i = 0;

susceptHeader = {'LightDur_h','Column','Period_min','GWS_Min360','GWS_Min60','AttenuationRatio','Label'};
maxSusRows = max(1, nCols * 250);
susceptTable = cell(maxSusRows, numel(susceptHeader));
suscept_i = 0;

% Track errors and write them out later
errors = cell(nCols, 2);
err_i  = 0;

% Upper bound for window-level rows
nW_max = max(1, ceil(durationHoursAll / stepHours) + 2);
P0_min_max = anchorBand_h(2) * 60;

sel_i  = struct();
full_i = struct();
qc_i   = struct();
win_i  = struct();

for mp = minPeriodsToRun
    key = sprintf('Min%d', mp);
    selReport.(key)    = cell(nCols, numel(selHeader));
    fullReport.(key)   = cell(nCols, numel(fullHeader));
    qcSelFull.(key)    = cell(nCols, numel(qcSelFullHeader));

    sel_i.(key)  = 0;
    full_i.(key) = 0;
    qc_i.(key)   = 0;

    Kmax = max(2, min(24, floor(P0_min_max / mp)));
    rowsPerWin = Kmax;
    maxWinRows = max(1, nCols * nW_max * rowsPerWin);

    selWinReport.(key) = cell(maxWinRows, numel(selWinHeader));
    win_i.(key) = 0;
end

%% ----------------------- OUTPUT ARRAYS (FULL LENGTH) ---------------------

resFull = struct(); remFull = struct();
resSel  = struct(); remSel  = struct();

for mp = minPeriodsToRun
    key = sprintf('Min%d', mp);
    resFull.(key) = NaN(N, numel(dataIdx));
    remFull.(key) = NaN(N, numel(dataIdx));
    resSel.(key)  = NaN(N, numel(dataIdx));
    remSel.(key)  = NaN(N, numel(dataIdx));
end

%% ----------------------- PROCESS EACH DATA COLUMN ------------------------

for c = 1:numel(dataIdx)
    colName = varNames{dataIdx(c)};
    fprintf('\nProcessing %d of %d: %s\n', c, numel(dataIdx), colName);

    try
        x0 = dataTable{:, dataIdx(c)};
        if ~isnumeric(x0)
            x0 = str2double(string(x0));
        end
        x0 = x0(:);

        % Light duration condition value for this dataset (usually constant)
        LightDur_h = light_value_repr(lightVec);

        % ------------------- Missingness and gap fill -------------------
        missingOriginal = isnan(x0) | ~isfinite(x0);
        missingFrac = sum(missingOriginal) / numel(x0);

        doInterpolate = (missingFrac <= missingFracThreshold) && (maxGapSamples > 0);

        xFill = x0;
        filledSamples = 0;
        if doInterpolate
            xFill = fillmissing(x0, 'pchip', 'MaxGap', maxGapSamples);
            filledSamples = sum(missingOriginal & ~isnan(xFill));
        end

        longestGapSamples = longest_nan_run(missingOriginal);
        [filledGapsCount, longGapsCount, unfilledSamples] = quantify_gaps(missingOriginal, isnan(xFill));

        gap_i = gap_i + 1;
        gapReport(gap_i, :) = {LightDur_h, colName, missingFrac, doInterpolate, filledSamples, ...
            filledGapsCount, longGapsCount, unfilledSamples, longestGapSamples};

        % ------------------- Baseline detrend -------------------
        base = movmedian(xFill, baselineWinSamples, 'omitnan', 'Endpoints', 'shrink');
        yDetrFull = xFill - base;

        validIdx = find(~isnan(yDetrFull) & isfinite(timeMinutesAll));
        if numel(validIdx) < minSamplesForAnalysis
            fprintf('  Too few valid samples for analysis. Outputs unchanged.\n');

            for mp = minPeriodsToRun
                key = sprintf('Min%d', mp);

                resFull.(key)(:, c) = x0;  remFull.(key)(:, c) = NaN(size(x0));
                resSel.(key)(:, c)  = x0;  remSel.(key)(:, c)  = NaN(size(x0));

                full_i.(key) = full_i.(key) + 1;
                fullReport.(key)(full_i.(key), :) = {LightDur_h, mp, colName, false, NaN, 0, 0, 'Too few valid samples'};

                sel_i.(key) = sel_i.(key) + 1;
                selReport.(key)(sel_i.(key), :)  = {LightDur_h, mp, colName, false, NaN, SUBTRACT_FUNDAMENTAL_IN_SELECTIVE, 0, '', 0, 0, 'Too few valid samples'};

                qc_i.(key) = qc_i.(key) + 1;
                qcSelFull.(key)(qc_i.(key), :)  = {LightDur_h, mp, colName, false, NaN, NaN, NaN, NaN, NaN, NaN, 'NA', 'Too few valid samples'};
            end

            anchor_i = anchor_i + 1;
            anchorReport(anchor_i, :) = {LightDur_h, colName, false, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'Too few valid samples'};
            continue;
        end

        t_min = timeMinutesAll(validIdx);
        t_hr  = t_min / 60;

        durCol_h = (max(t_min, [], 'omitnan') - min(t_min, [], 'omitnan')) / 60;

        y = yDetrFull(validIdx);
        y = y - median(y, 'omitnan');

        % ------------------- STEP 1: Regression anchor scan 22-28 h -------------------

        periods_h = (anchorBand_h(1):anchorStep_h:anchorBand_h(2))';
        [P0_h, maxDeltaR2, peakZ, A0, ampSNR, deltaCurve] = regression_anchor_scan(t_hr, y, periods_h);
        cyclesAtP0 = durCol_h / max(P0_h, eps);

        % ------------------- STEP 2: Block-shuffle surrogate p-value -------------------

        if numel(y) < 4*blockLenSamples
            blk = max(4, floor(numel(y)/4));
        else
            blk = blockLenSamples;
        end

        surrMax = NaN(Nsurr, 1);
        for s = 1:Nsurr
            yS = block_shuffle_surrogate(y, blk);
            surrMax(s) = regression_anchor_maxDeltaR2(t_hr, yS, periods_h);
        end
        pBlock = (1 + sum(surrMax >= maxDeltaR2)) / (1 + Nsurr);

        % ------------------- AnchorOK gate -------------------

        reasons = {};
        warnings = {};

        finiteOK = isfinite(P0_h) && isfinite(maxDeltaR2) && isfinite(pBlock) && isfinite(cyclesAtP0);
        if ~finiteOK
            reasons = [reasons, {'Non-finite anchor statistics'}];
        end

        sigOK = (pBlock < alphaAnchor);
        if ~sigOK
            reasons = [reasons, {sprintf('Block-shuffle p too large (%.3g >= %.3g)', pBlock, alphaAnchor)}];
        end

        cyclesOK = (cyclesAtP0 >= minCyclesForAnchor);
        if ~cyclesOK
            reasons = [reasons, {sprintf('Too few cycles (%.2f < %.1f)', cyclesAtP0, minCyclesForAnchor)}];
        end

        deltaOK = (maxDeltaR2 >= minDeltaR2);
        if ~deltaOK
            reasons = [reasons, {sprintf('DeltaR2 too small (%.4g < %.4g)', maxDeltaR2, minDeltaR2)}];
        end

        edgeOK = true;
        if useEdgeGuard
            edgeOK = (P0_h > (anchorBand_h(1) + edgeMargin_h)) && (P0_h < (anchorBand_h(2) - edgeMargin_h));
            if ~edgeOK
                reasons = [reasons, {sprintf('Period too close to band edge (Period=%.3f h)', P0_h)}];
            end
        end

        if isfinite(peakZ) && peakZ < PeakZ_warn
            warnings = [warnings, {sprintf('QC: broad peak (PeakZ=%.2f < %.2f)', peakZ, PeakZ_warn)}]; %#ok<*AGROW>
        end
        if isfinite(ampSNR) && ampSNR < SNR_warn
            warnings = [warnings, {sprintf('QC: low SNR (AmpSNR=%.2f < %.2f)', ampSNR, SNR_warn)}];
        end

        anchorOK = finiteOK && sigOK && cyclesOK && deltaOK && edgeOK;

        if isempty(reasons)
            noteAnchor = 'Anchor accepted (pBlock + minimal viability).';
        else
            noteAnchor = strjoin(reasons, '; ');
        end
        if ~isempty(warnings)
            noteAnchor = [noteAnchor ' | ' strjoin(warnings, '; ')];
        end

        fprintf('  Anchor: %s | Period=%.3f h | maxDeltaR2=%.4g | PeakZ=%.2f | Amp=%.4g | SNR=%.2f | cycles=%.2f | pBlock=%.3g\n', ...
            tf(anchorOK), P0_h, maxDeltaR2, peakZ, A0, ampSNR, cyclesAtP0, pBlock);

        anchor_i = anchor_i + 1;
        anchorReport(anchor_i, :) = {LightDur_h, colName, anchorOK, P0_h, maxDeltaR2, peakZ, A0, ampSNR, durCol_h, cyclesAtP0, pBlock, noteAnchor};

        % ------------------- Anchor figure -------------------

        if SAVE_ANCHOR_FIGS
            try
                fig = figure('Visible','off','Color','w','Position',[100 100 980 420]);
                tiledlayout(fig, 1, 2, 'Padding','compact','TileSpacing','compact');

                nexttile;
                plot(periods_h, deltaCurve, 'LineWidth', 1.2);
                hold on; xline(P0_h, 'LineWidth', 1.2); hold off;
                set(gca, 'Box','off', 'FontName','Times New Roman', 'TickDir','out', ...
                    'XGrid','off','YGrid','off');
                xlabel('Period (h)', 'FontWeight','bold', 'FontName','Times New Roman');
                ylabel('DeltaR^2',    'FontWeight','bold', 'FontName','Times New Roman');
                title(sprintf('%s | Period %.2f h', colName, P0_h), 'FontName','Times New Roman');

                nexttile;
                histogram(surrMax, 25);
                hold on; xline(maxDeltaR2, 'LineWidth', 1.2); hold off;
                set(gca, 'Box','off', 'FontName','Times New Roman', 'TickDir','out', ...
                    'XGrid','off','YGrid','off');
                xlabel('Surrogate max DeltaR^2 (22-28 h)', 'FontWeight','bold', 'FontName','Times New Roman');
                ylabel('Count', 'FontWeight','bold', 'FontName','Times New Roman');
                title(sprintf('Block-shuffle null | p=%.3g', pBlock), 'FontName','Times New Roman');

                outFig = fullfile(figAnchorFolder, sprintf('Anchor_%s.jpg', sanitise_filename(colName)));
                print(fig, outFig, '-djpeg', '-r600');
                close(fig);
            catch
            end
        end

        % ------------------- If anchor rejected: pass through -------------------
        
        if ~anchorOK
            for mp = minPeriodsToRun
                key = sprintf('Min%d', mp);

                resFull.(key)(:, c) = x0;  remFull.(key)(:, c) = NaN(size(x0));
                resSel.(key)(:, c)  = x0;  remSel.(key)(:, c)  = NaN(size(x0));

                full_i.(key) = full_i.(key) + 1;
                fullReport.(key)(full_i.(key), :) = {LightDur_h, mp, colName, false, P0_h, 0, 0, 'Anchor rejected (no subtraction)'};

                sel_i.(key) = sel_i.(key) + 1;
                selReport.(key)(sel_i.(key), :)  = {LightDur_h, mp, colName, false, P0_h, SUBTRACT_FUNDAMENTAL_IN_SELECTIVE, 0, '', 0, 0, 'Anchor rejected (no subtraction)'};

                qc_i.(key) = qc_i.(key) + 1;
                qcSelFull.(key)(qc_i.(key), :)  = {LightDur_h, mp, colName, false, P0_h, NaN, NaN, NaN, NaN, NaN, 'NA', 'Anchor rejected'};
            end
            continue;
        end

        % ------------------- Windowing for coupling (dynamic) -------------------

        winHours = min(max(72, round(durCol_h/2)), 168);

        tMinAll = min(t_min); tMaxAll = max(t_min);
        [winStarts_min, winEnds_min] = make_windows(tMinAll, tMaxAll, winHours*60, stepHours*60);
        nW = numel(winStarts_min);

        fprintf('  Coupling windows: %g h | step %g h | nW=%d\n', winHours, stepHours, nW);

        % ------------------- Fundamental per-window fit at P0 -------------------

        A1_w   = NaN(nW,1);
        Phi1_w = NaN(nW,1);
        P0_min = P0_h * 60;

        for w = 1:nW
            wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
            if sum(wMask) < minSamplesPerWindow
                continue;
            end
            tmw = t_min(wMask);
            yw  = y(wMask);

            [A1, phi1, R2w] = fit_sinusoid_amp_phase_R2_with_linear(tmw, yw, P0_min);

            A1_w(w)   = A1;
            Phi1_w(w) = phi1;

            for mp = minPeriodsToRun
                key = sprintf('Min%d', mp);
                win_i.(key) = win_i.(key) + 1;
                selWinReport.(key)(win_i.(key), :) = {LightDur_h, mp, colName, w, winStarts_min(w)/60, winEnds_min(w)/60, P0_h, 1, A1, phi1, R2w, P0_h, 1, 'Fundamental (fixed Period, +linear)'};
            end
        end

        okFund = isfinite(A1_w) & isfinite(Phi1_w);
        nFund  = sum(okFund);

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

        % Prepare scalogram filterbank (common) if enabled
        if SAVE_SCALOGRAMS
            if ~exist(figScaloFolder, 'dir'); mkdir(figScaloFolder); end
            if ~exist(figScaloResFolder, 'dir'); mkdir(figScaloResFolder); end
            if ~exist(figScaloRemFolder, 'dir'); mkdir(figScaloRemFolder); end

            scaloFB = cwtfilterbank('SignalLength', N, ...
                'SamplingPeriod', minutes(TsMinutes), ...
                'PeriodLimits', [minutes(SCALO_MIN_PERIOD_MIN), minutes(SCALO_MAX_PERIOD_MIN)], ...
                'Wavelet', 'amor');
        end

        % ------------------- Per-minPeriod passes -------------------

        tmpFull = struct();
        tmpVarExplFull = struct();
        tmpK = struct();

        for mp = minPeriodsToRun
            key = sprintf('Min%d', mp);
            minPeriodInterestMin = mp;

            % ------------------- Harmonic ladder depth -------------------

            K = floor(P0_min / minPeriodInterestMin);
            K = max(2, min(K, 24)); % cap
            kList = 2:K;

            % ------------------- Harmonic candidates (per window) -------------------

            Ak_w     = NaN(nW, numel(kList));
            Phik_w   = NaN(nW, numel(kList));
            ratio_w  = NaN(nW, numel(kList));
            Pbest_h  = NaN(nW, numel(kList));
            R2k_w    = NaN(nW, numel(kList));

            for w = 1:nW
                wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
                if sum(wMask) < minSamplesPerWindow || ~isfinite(A1_w(w)) || ~isfinite(Phi1_w(w))
                    continue;
                end
                tmw = t_min(wMask);
                yw  = y(wMask);

                for j = 1:numel(kList)
                    k = kList(j);

                    targetP_min = P0_min / k;
                    scanP_min = linspace(targetP_min*(1-peakSearchFrac), targetP_min*(1+peakSearchFrac), 19);

                    bestR2  = -Inf;
                    bestA   = NaN;
                    bestPhi = NaN;
                    bestP   = NaN;

                    for sp = 1:numel(scanP_min)
                        [Akk, phikk, R2kk] = fit_sinusoid_amp_phase_R2_with_linear(tmw, yw, scanP_min(sp));
                        if isfinite(R2kk) && R2kk > bestR2
                            bestR2  = R2kk;
                            bestA   = Akk;
                            bestPhi = phikk;
                            bestP   = scanP_min(sp);
                        end
                    end

                    if ~isfinite(bestR2)
                        continue;
                    end

                    Ak_w(w, j)    = bestA;
                    Phik_w(w, j)  = bestPhi;
                    ratio_w(w, j) = bestP / targetP_min;
                    Pbest_h(w, j) = bestP / 60;
                    R2k_w(w, j)   = bestR2;

                    win_i.(key) = win_i.(key) + 1;
                    selWinReport.(key)(win_i.(key), :) = {LightDur_h, mp, colName, w, winStarts_min(w)/60, winEnds_min(w)/60, ...
                        P0_h, k, bestA, bestPhi, bestR2, bestP/60, ratio_w(w, j), 'Harmonic candidate (scan fit, +linear)'};
                end
            end

            % ------------------- Harmonic classification (Selective) -------------------

            harmonicLikely = false(1, numel(kList));

            for j = 1:numel(kList)
                k = kList(j);

                okW = isfinite(A1_w) & isfinite(Ak_w(:, j)) & isfinite(Phi1_w) & isfinite(Phik_w(:, j)) & isfinite(ratio_w(:, j));
                nOk = sum(okW);

                if nOk < minWindowsForDecision
                    continue;
                end

                ratios = ratio_w(okW, j);
                ratioMad = mad(ratios - 1, 1);
                ratioPass = ratioMad <= ratioTol;

                del = wrapToPi_local(Phik_w(okW, j) - k * Phi1_w(okW));
                Rplv = abs(mean(exp(1i * del)));
                Rpass = isfinite(Rplv) && (Rplv >= Rcrit);

                if strcmp(regime, 'long')
                    rho = corr(A1_w(okW), Ak_w(okW, j), 'Type', 'Spearman', 'Rows', 'complete');
                    rhoPass = isfinite(rho) && (abs(rho) >= rhoCrit);
                    harmonicLikely(j) = (ratioPass + rhoPass + Rpass) >= 2;
                elseif strcmp(regime, 'short')
                    A1_med = median(A1_w(okW), 'omitnan');
                    Ak_med = median(Ak_w(okW, j), 'omitnan');
                    ampPass = isfinite(A1_med) && A1_med > 0 && isfinite(Ak_med) && (Ak_med >= ampFracCritShort * A1_med);
                    harmonicLikely(j) = ratioPass && Rpass && ampPass;
                end
            end

            kLikely = kList(harmonicLikely);

            % ------------------- Selective subtraction -------------------

            if SUBTRACT_FUNDAMENTAL_IN_SELECTIVE
                kSubSel = unique([1, kLikely(:)'], 'stable');
            else
                kSubSel = kLikely(:)'; %#ok<*UNRCH>
            end

            yFit = y;

            if isempty(kSubSel)
                removedSel = zeros(size(yFit));
                residualSel = yFit;
                varExplSel = 0;
                noteSel = sprintf('Anchor OK (Period=%.3f h) but no harmonics classified (%s regime, nWin=%d).', P0_h, regime, nFund);
            else
                [X0, Xh] = build_drift_and_harmonics_design(t_hr, t_min, P0_min, kSubSel);
                X = [X0, Xh];
                beta = X \ yFit;
                yhatDrift = X0 * beta(1:size(X0,2));
                yhatTotal = X  * beta;

                removedSel  = yhatTotal - yhatDrift;
                residualSel = yFit - removedSel;

                varExplSel = variance_explained(yFit, residualSel);
                noteSel = sprintf('Selective removed k=%s (%s regime, nWin=%d).', mat2str(kSubSel), regime, nFund);
            end

            [rSelFull, remSelFull_] = map_back_to_full(x0, validIdx, residualSel, removedSel);
            resSel.(key)(:, c) = rSelFull;
            remSel.(key)(:, c) = remSelFull_;

            sel_i.(key) = sel_i.(key) + 1;
            selReport.(key)(sel_i.(key), :) = {LightDur_h, mp, colName, true, P0_h, SUBTRACT_FUNDAMENTAL_IN_SELECTIVE, ...
                numel(kSubSel), mat2str(kSubSel), varExplSel, nFund, noteSel};

            % ------------------- FullLadder subtraction -------------------

            kSubFull = 1:K;

            [X0F, XhF] = build_drift_and_harmonics_design(t_hr, t_min, P0_min, kSubFull);
            XF = [X0F, XhF];
            betaF = XF \ yFit;
            yhatDriftF = X0F * betaF(1:size(X0F,2));
            yhatTotalF = XF  * betaF;

            removedFull  = yhatTotalF - yhatDriftF;
            residualFull = yFit - removedFull;
            varExplFullV = variance_explained(yFit, residualFull);

            [rFullFull, remFullFull_] = map_back_to_full(x0, validIdx, residualFull, removedFull);
            resFull.(key)(:, c) = rFullFull;
            remFull.(key)(:, c) = remFullFull_;

            full_i.(key) = full_i.(key) + 1;
            fullReport.(key)(full_i.(key), :) = {LightDur_h, mp, colName, true, P0_h, K, varExplFullV, ...
                sprintf('FullLadder removed k=1..%d (minPeriod=%g min)', K, minPeriodInterestMin)};

            fprintf('  %s | Selective VarExpl=%.3f | FullLadder VarExpl=%.3f (K=%d)\n', key, varExplSel, varExplFullV, K);

            % ------------------- Selective vs FullLadder QC (same minPeriod) -------------------

            rSel_v  = residualSel(:);
            rFull_v = residualFull(:);

            rmsDiff = sqrt(mean((rSel_v - rFull_v).^2, 'omitnan'));
            denom   = sqrt(mean((yFit).^2, 'omitnan'));
            nRMS    = rmsDiff / (denom + eps);

            dVar = varExplFullV - varExplSel;

            risk = 'Low';
            if (nRMS >= nRMS_high) || (dVar >= dVar_high)
                risk = 'High';
            elseif (nRMS >= nRMS_low) || (dVar >= dVar_low)
                risk = 'Moderate';
            end

            qc_i.(key) = qc_i.(key) + 1;
            qcSelFull.(key)(qc_i.(key), :) = {LightDur_h, mp, colName, true, P0_h, ...
                varExplSel, varExplFullV, dVar, rmsDiff, nRMS, risk, sprintf('Regime=%s, nWin=%d, K=%d', regime, nFund, K)};

            % ------------------- Scalograms (FullLadder + Selective) -------------------

            if SAVE_SCALOGRAMS
                try
                    % FullLadder residual/removed
                    sResF = rFullFull(:);    sResF(~isfinite(sResF)) = 0;
                    sRemF = remFullFull_(:); sRemF(~isfinite(sRemF)) = 0;

                    [wtResF, periodsResF] = cwt(sResF, 'FilterBank', scaloFB);
                    [wtRemF, periodsRemF] = cwt(sRemF, 'FilterBank', scaloFB);

                    pH_resF = hours(periodsResF);
                    pH_remF = hours(periodsRemF);

                    safeCol = sanitise_filename(colName);

                    % Residual scalogram (FullLadder)
                    figS1 = figure('Visible','off','Color','w','Position',[100 100 980 640]);
                    pcolor(time_day, pH_resF, abs(wtResF));
                    shading interp;
                    colormap jet;
                    colorbar;
                    caxis auto;
                    set(gca,'box','off','YTick',SCALO_YTICKS_H,'TickDir','out','FontName','Times New Roman', ...
                        'XGrid','off','YGrid','off');
                    hold on;
                    for kLine = 1:numel(condChangeIdx)
                        boundaryRow = condChangeIdx(kLine) + 1;
                        if boundaryRow >= 1 && boundaryRow <= numel(time_day)
                            x_line = time_day(boundaryRow);
                            plot([x_line x_line], [min(pH_resF) max(pH_resF)], 'w:', 'LineWidth', 1.5);
                        end
                    end
                    hold off;
                    xlabel('Time (days)', 'FontWeight','bold', 'FontName','Times New Roman');
                    ylabel('Period (hr)', 'FontWeight','bold', 'FontName','Times New Roman');
                    title(sprintf('Scalogram - Residual | FullLadder %s | %s', key, colName), 'FontName','Times New Roman');
                    outResF = fullfile(figScaloResFolder, sprintf('Scalogram_Res_FullLadder_%s_%s.jpg', key, safeCol));
                    print(figS1, outResF, '-djpeg', '-r600');
                    close(figS1);

                    % Removed scalogram (FullLadder)
                    figS2 = figure('Visible','off','Color','w','Position',[100 100 980 640]);
                    pcolor(time_day, pH_remF, abs(wtRemF));
                    shading interp;
                    colormap jet;
                    colorbar;
                    caxis auto;
                    set(gca,'box','off','YTick',SCALO_YTICKS_H,'TickDir','out','FontName','Times New Roman', ...
                        'XGrid','off','YGrid','off');
                    hold on;
                    for kLine = 1:numel(condChangeIdx)
                        boundaryRow = condChangeIdx(kLine) + 1;
                        if boundaryRow >= 1 && boundaryRow <= numel(time_day)
                            x_line = time_day(boundaryRow);
                            plot([x_line x_line], [min(pH_remF) max(pH_remF)], 'w:', 'LineWidth', 1.5);
                        end
                    end
                    hold off;
                    xlabel('Time (days)', 'FontWeight','bold', 'FontName','Times New Roman');
                    ylabel('Period (hr)', 'FontWeight','bold', 'FontName','Times New Roman');
                    title(sprintf('Scalogram - Removed | FullLadder %s | %s', key, colName), 'FontName','Times New Roman');
                    outRemF = fullfile(figScaloRemFolder, sprintf('Scalogram_Rem_FullLadder_%s_%s.jpg', key, safeCol));
                    print(figS2, outRemF, '-djpeg', '-r600');
                    close(figS2);

                    % Selective scalograms (mirroring FullLadder)
                    if SAVE_SCALOGRAMS_SELECTIVE
                        sResS = rSelFull(:);      sResS(~isfinite(sResS)) = 0;
                        sRemS = remSelFull_(:);   sRemS(~isfinite(sRemS)) = 0;

                        [wtResS, periodsResS] = cwt(sResS, 'FilterBank', scaloFB);
                        [wtRemS, periodsRemS] = cwt(sRemS, 'FilterBank', scaloFB);

                        pH_resS = hours(periodsResS);
                        pH_remS = hours(periodsRemS);

                        % Residual scalogram (Selective)
                        figSS1 = figure('Visible','off','Color','w','Position',[100 100 980 640]);
                        pcolor(time_day, pH_resS, abs(wtResS));
                        shading interp;
                        colormap jet;
                        colorbar;
                        caxis auto;
                        set(gca,'box','off','YTick',SCALO_YTICKS_H,'TickDir','out','FontName','Times New Roman', ...
                            'XGrid','off','YGrid','off');
                        hold on;
                        for kLine = 1:numel(condChangeIdx)
                            boundaryRow = condChangeIdx(kLine) + 1;
                            if boundaryRow >= 1 && boundaryRow <= numel(time_day)
                                x_line = time_day(boundaryRow);
                                plot([x_line x_line], [min(pH_resS) max(pH_resS)], 'w:', 'LineWidth', 1.5);
                            end
                        end
                        hold off;
                        xlabel('Time (days)', 'FontWeight','bold', 'FontName','Times New Roman');
                        ylabel('Period (hr)', 'FontWeight','bold', 'FontName','Times New Roman');
                        title(sprintf('Scalogram - Residual | Selective %s | %s', key, colName), 'FontName','Times New Roman');
                        outResS = fullfile(figScaloResFolder, sprintf('Scalogram_Res_Selective_%s_%s.jpg', key, safeCol));
                        print(figSS1, outResS, '-djpeg', '-r600');
                        close(figSS1);

                        % Removed scalogram (Selective)
                        figSS2 = figure('Visible','off','Color','w','Position',[100 100 980 640]);
                        pcolor(time_day, pH_remS, abs(wtRemS));
                        shading interp;
                        colormap jet;
                        colorbar;
                        caxis auto;
                        set(gca,'box','off','YTick',SCALO_YTICKS_H,'TickDir','out','FontName','Times New Roman', ...
                            'XGrid','off','YGrid','off');
                        hold on;
                        for kLine = 1:numel(condChangeIdx)
                            boundaryRow = condChangeIdx(kLine) + 1;
                            if boundaryRow >= 1 && boundaryRow <= numel(time_day)
                                x_line = time_day(boundaryRow);
                                plot([x_line x_line], [min(pH_remS) max(pH_remS)], 'w:', 'LineWidth', 1.5);
                            end
                        end
                        hold off;
                        xlabel('Time (days)', 'FontWeight','bold', 'FontName','Times New Roman');
                        ylabel('Period (hr)', 'FontWeight','bold', 'FontName','Times New Roman');
                        title(sprintf('Scalogram - Removed | Selective %s | %s', key, colName), 'FontName','Times New Roman');
                        outRemS = fullfile(figScaloRemFolder, sprintf('Scalogram_Rem_Selective_%s_%s.jpg', key, safeCol));
                        print(figSS2, outRemS, '-djpeg', '-r600');
                        close(figSS2);
                    end
                catch
                end
            end

            % Store for cross-min comparisons and susceptibility (FullLadder)
            tmpFull.(key) = residualFull(:);
            tmpVarExplFull.(key) = varExplFullV;
            tmpK.(key) = K;
        end

        % ------------------- Cross-min QC (FullLadder Min360 vs Min60) -------------------

        if isfield(tmpFull, 'Min360') && isfield(tmpFull, 'Min60')
            r360 = tmpFull.Min360;
            r60  = tmpFull.Min60;

            rmsCross = sqrt(mean((r60 - r360).^2, 'omitnan'));
            denomCross = sqrt(mean((y).^2, 'omitnan'));
            nRMS_cross = rmsCross / (denomCross + eps);

            dVar_cross = tmpVarExplFull.Min60 - tmpVarExplFull.Min360;

            rec = 'Min360';
            if (nRMS_cross >= nRMS_high) || (dVar_cross >= dVar_high)
                rec = 'Min60';
            elseif (nRMS_cross >= nRMS_low) || (dVar_cross >= dVar_low)
                rec = 'Review';
            end

            cross_i = cross_i + 1;
            crossMinQC(cross_i, :) = {LightDur_h, colName, true, P0_h, ...
                tmpK.Min360, tmpK.Min60, tmpVarExplFull.Min360, tmpVarExplFull.Min60, dVar_cross, rmsCross, nRMS_cross, rec, ...
                sprintf('Thresholds: nRMS(%.2f/%.2f), dVar(%.2f/%.2f)', nRMS_low, nRMS_high, dVar_low, dVar_high)};
        end

        % ------------------- Harmonic susceptibility periods (FullLadder) -------------------

        if SAVE_SCALOGRAMS && exist('scaloFB', 'var') && isfield(tmpFull, 'Min360') && isfield(tmpFull, 'Min60')
            try
                s360 = resFull.Min360(:, c); s360(~isfinite(s360)) = 0;
                s60  = resFull.Min60(:,  c); s60(~isfinite(s60))  = 0;

                [wt360, per360] = cwt(s360, 'FilterBank', scaloFB);
                [wt60,  ~    ]  = cwt(s60,  'FilterBank', scaloFB);

                perMin = minutes(per360);

                gws360 = mean(abs(wt360).^2, 2, 'omitnan');
                gws60  = mean(abs(wt60 ).^2, 2, 'omitnan');

                for iP = 1:numel(perMin)
                    ratio = gws60(iP) / (gws360(iP) + eps);
                    label = 'Ambiguous';
                    if ratio < 0.30
                        label = 'Harmonic-susceptible';
                    elseif ratio > 0.70
                        label = 'Harmonic-robust';
                    end

                    suscept_i = suscept_i + 1;
                    if suscept_i > size(susceptTable, 1)
                        susceptTable = [susceptTable; cell(maxSusRows, numel(susceptHeader))];
                    end
                    susceptTable(suscept_i, :) = {LightDur_h, colName, perMin(iP), gws360(iP), gws60(iP), ratio, label};
                end
            catch
            end
        end

    catch ME
        fprintf('  Error: %s\n', ME.message);
        err_i = err_i + 1;
        errors(err_i, :) = {fileName, sprintf('Column %s: %s', colName, ME.message)};
    end
end

%% ------------------------ TRIM PREALLOCATED LOGS --------------------------

gapReport    = gapReport(1:gap_i, :);
anchorReport = anchorReport(1:anchor_i, :);
crossMinQC   = crossMinQC(1:cross_i, :);
susceptTable = susceptTable(1:suscept_i, :);
errors       = errors(1:err_i, :);

for mp = minPeriodsToRun
    key = sprintf('Min%d', mp);
    selReport.(key)    = selReport.(key)(1:sel_i.(key), :);
    fullReport.(key)   = fullReport.(key)(1:full_i.(key), :);
    qcSelFull.(key)    = qcSelFull.(key)(1:qc_i.(key), :);
    selWinReport.(key) = selWinReport.(key)(1:win_i.(key), :);
end

%% -------------------- WRITE TIME-SERIES EXCEL FILES ----------------------

try
    timeOut  = dataTable{:, timeIdx};
    lightOut = dataTable{:, lightIdx};

    for mp = minPeriodsToRun
        key = sprintf('Min%d', mp);

        write_output_excel(tsResidualFolder, fileName, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            resFull.(key), sprintf('Residual_FullLadder_%s', key));

        write_output_excel(tsRemovedFolder, fileName, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            remFull.(key), sprintf('Removed_FullLadder_%s', key));

        write_output_excel(tsResidualFolder, fileName, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            resSel.(key), sprintf('Residual_Selective_%s', key));

        write_output_excel(tsRemovedFolder, fileName, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            remSel.(key), sprintf('Removed_Selective_%s', key));
    end
catch ME
    fprintf('Output write failed: %s\n', ME.message);
end

%% ------------------------ WRITE REPORT WORKBOOKS -------------------------

gapTable      = cell2table(gapReport,      'VariableNames', gapReportHeader);
anchorTable   = cell2table(anchorReport,   'VariableNames', anchorHeader);
excludedTable = cell2table(excludedReport, 'VariableNames', excludedHeader);
errorsTable   = cell2table(errors,         'VariableNames', {'File','Error'});

T_sel  = struct(); T_full = struct(); T_qc = struct(); T_win = struct();
for mp = minPeriodsToRun
    key = sprintf('Min%d', mp);
    T_sel.(key)  = cell2table(selReport.(key),    'VariableNames', selHeader);
    T_full.(key) = cell2table(fullReport.(key),   'VariableNames', fullHeader);
    T_qc.(key)   = cell2table(qcSelFull.(key),    'VariableNames', qcSelFullHeader);
    T_win.(key)  = cell2table(selWinReport.(key), 'VariableNames', selWinHeader);
end

T_cross = cell2table(crossMinQC,    'VariableNames', crossMinHeader);
T_sus   = cell2table(susceptTable,  'VariableNames', susceptHeader);

%% ------------------------ BUILD RECOMMENDATION SHEET ----------------------

recHeader = {'LightDur_h','Column','AnchorOK','Period_hours','RecommendedResidual','Workbook','SignalColumnHeader','Reason'};
recRows = cell(0, numel(recHeader));

for c = 1:numel(dataIdx)
    colName = varNames{dataIdx(c)};

    aMask = strcmp(string(anchorTable.Column), string(colName));
    if any(aMask)
        aRow = anchorTable(find(aMask, 1, 'first'), :);
        anchorOK = logical(aRow.AnchorOK);
        P0_h = aRow.Period_hours;
        LightDur_h = aRow.LightDur_h;
    else
        anchorOK = false;
        P0_h = NaN;
        LightDur_h = NaN;
    end

    chosenMinKey = 'Min360';
    crossReason = ''; %#ok<*NASGU>
    crossMetrics = '';

    if anchorOK
        cMask = strcmp(string(T_cross.Column), string(colName));
        if any(cMask)
            cRow = T_cross(find(cMask, 1, 'first'), :);
            nRMS_cross = cRow.nRMS_60vs360;
            dVar_cross = cRow.DeltaVarExpl_60minus360;

            if (nRMS_cross >= nRMS_high) || (dVar_cross >= dVar_high)
                chosenMinKey = 'Min60';
                crossReason = sprintf('Chose Min60 because cross-min nRMS_60vs360=%.3f (>=%.2f) and/or DeltaVarExpl_60minus360=%.3f (>=%.2f).', ...
                    nRMS_cross, nRMS_high, dVar_cross, dVar_high);
            elseif (nRMS_cross >= nRMS_low) || (dVar_cross >= dVar_low)
                chosenMinKey = 'Min360';
                crossReason = sprintf('Chose Min360 (borderline cross-min): nRMS_60vs360=%.3f (>=%.2f) or DeltaVarExpl_60minus360=%.3f (>=%.2f), but not exceeding high thresholds.', ...
                    nRMS_cross, nRMS_low, dVar_cross, dVar_low);
            else
                chosenMinKey = 'Min360';
                crossReason = sprintf('Chose Min360 because cross-min differences are small: nRMS_60vs360=%.3f (<%.2f) and DeltaVarExpl_60minus360=%.3f (<%.2f).', ...
                    nRMS_cross, nRMS_low, dVar_cross, dVar_low);
            end
            crossMetrics = sprintf('CrossMin: nRMS_60vs360=%.3f, DeltaVarExpl_60minus360=%.3f', nRMS_cross, dVar_cross);
        else
            chosenMinKey = 'Min360';
            crossReason = 'Chose Min360 (no cross-min QC row available).';
            crossMetrics = 'CrossMin: NA';
        end
    else
        crossReason = 'Anchor rejected: no subtraction performed. All residuals are pass-through for this mouse.';
        crossMetrics = 'CrossMin: NA';
    end

    mode = 'Selective';
    withinReason = '';
    withinMetrics = '';

    if anchorOK
        qcTbl = T_qc.(chosenMinKey);
        qMask = strcmp(string(qcTbl.Column), string(colName));
        if any(qMask)
            qRow = qcTbl(find(qMask, 1, 'first'), :);
            nRMS = qRow.nRMS;
            dVar = qRow.DeltaVarExpl_FullMinusSel;

            if (nRMS < nRMS_low) && (dVar < dVar_low)
                mode = 'Selective';
                withinReason = sprintf('Chose Selective because within-min disagreement is low: nRMS=%.3f (<%.2f) and DeltaVarExpl_FullMinusSel=%.3f (<%.2f).', ...
                    nRMS, nRMS_low, dVar, dVar_low);
            elseif (nRMS >= nRMS_high) || (dVar >= dVar_high)
                mode = 'FullLadder';
                withinReason = sprintf('Chose FullLadder because Selective vs FullLadder disagreement is high: nRMS=%.3f (>=%.2f) and/or DeltaVarExpl_FullMinusSel=%.3f (>=%.2f).', ...
                    nRMS, nRMS_high, dVar, dVar_high);
            else
                mode = 'Selective';
                withinReason = sprintf('Chose Selective (default) because disagreement is moderate: nRMS=%.3f and DeltaVarExpl_FullMinusSel=%.3f (not exceeding high thresholds).', ...
                    nRMS, dVar);
            end
            withinMetrics = sprintf('WithinMin(%s): nRMS=%.3f, DeltaVarExpl_FullMinusSel=%.3f', chosenMinKey, nRMS, dVar);
        else
            mode = 'Selective';
            withinReason = sprintf('Chose Selective (no within-min QC row available for %s).', chosenMinKey);
            withinMetrics = sprintf('WithinMin(%s): NA', chosenMinKey);
        end
    else
        mode = 'Selective';
        withinReason = 'Anchor rejected: defaulting to Selective_Min360 label (pass-through).';
        withinMetrics = 'WithinMin: NA';
    end

    recResidual = sprintf('%s_%s', mode, chosenMinKey);

    if strcmp(mode, 'FullLadder')
        prefix = sprintf('Residual_FullLadder_%s', chosenMinKey);
    else
        prefix = sprintf('Residual_Selective_%s', chosenMinKey);
    end

    workbook = fullfile('TimeSeries', 'Residual', sprintf('%s_%s', prefix, fileName));
    sigHeader = sprintf('%s_%s', prefix, colName);

    reason = strjoin({crossReason, withinReason, crossMetrics, withinMetrics}, ' ');

    recRows(end+1,:) = {LightDur_h, colName, anchorOK, P0_h, recResidual, workbook, sigHeader, reason}; %#ok<*SAGROW>
end

recommendationTable = cell2table(recRows, 'VariableNames', recHeader);

%% ---------------------------- README + INDEX -----------------------------

readmeLines = {
    'Harmonic removal reports (Summary workbook)'
    ''
    'Overview:'
    '- Anchor detection scans 22-28 h and uses block-shuffle p-value gating.'
    '- FullLadder removes k=1..K down to Min360 and Min60.'
    '- Selective removes k that behave like harmonics across windows (plus k=1).'
    '- Time-series outputs are in TimeSeries/Residual and TimeSeries/Removed.'
    '- Scalogram outputs are in Figures_Scalograms/Residual and Figures_Scalograms/Removed.'
    ''
    'Key QC interpretation:'
    'Anchor_Report:'
    '- AnchorOK=1 means supported circadian-scale component (pBlock<0.05 plus viability guards).'
    '- MaxDeltaR2 higher means stronger circadian-like sinusoid beyond drift.'
    '- PeakZ and AmpSNR are QC flags only (low means broad/noisy anchor).'
    ''
    'Selective_vs_FullLadder_QC (within Min):'
    sprintf('- nRMS < %.2f low disagreement; %.2f-%.2f moderate; >= %.2f high.', nRMS_low, nRMS_low, nRMS_high, nRMS_high)
    sprintf('- DeltaVarExpl_FullMinusSel < %.2f small; %.2f-%.2f moderate; >= %.2f large.', dVar_low, dVar_low, dVar_high, dVar_high)
    ''
    'CrossMinPeriod_QC (FullLadder Min60 vs Min360):'
    '- Large nRMS_60vs360 and DeltaVarExpl_60minus360 indicates high-order harmonic contribution.'
    ''
    'Recommendation sheet:'
    '- Picks one of the 4 residuals per mouse using CrossMinPeriod_QC (to choose Min360 vs Min60) and'
    '  Selective_vs_FullLadder_QC (to choose Selective vs FullLadder) with the thresholds above.'
    ''
    'Excluded_Columns sheet:'
    '- Records any Data columns excluded via the optional exclude selection step.'
    ''
    'Harmonic_Susceptibility_Periods (Detail workbook):'
    '- AttenuationRatio = GWS_Min60 / GWS_Min360.'
    '- <0.30 harmonic-susceptible; 0.30-0.70 ambiguous; >0.70 harmonic-robust.'
    ''
    'Errors sheet:'
    '- Captures any caught per-mouse exceptions during processing.'
    ''
    'Figures:'
    '- Anchor figures are saved as JPG at 600 DPI in Figures_Anchor/.'
    '- Scalograms (if enabled) are saved as JPG at 600 DPI in Figures_Scalograms/.'};

indexRows = {
    'Sheet','Contents'
    'README','How to interpret outputs and QC'
    'Index','Sheet map'
    'Excluded_Columns','Columns excluded by user at selection stage'
    'GapFill_Report','Missingness and interpolation actions per mouse'
    'Anchor_Report','Anchor period detection results and acceptance gating'
    'FullLadder_Min360','FullLadder removal summary at Min360'
    'FullLadder_Min60','FullLadder removal summary at Min60'
    'Selective_Min360','Selective removal summary at Min360'
    'Selective_Min60','Selective removal summary at Min60'
    'QC_SelVsFull_Min360','Selective vs FullLadder QC at Min360'
    'QC_SelVsFull_Min60','Selective vs FullLadder QC at Min60'
    'QC_CrossMin_360vs60','FullLadder cross-min QC'
    'Recommendation','Per-mouse recommended residual (one of the 4), with reason/metrics'
    'Errors','Any per-mouse processing exceptions caught during the run'};

%% ------------------------ WRITE SUMMARY WORKBOOK -------------------------

writecell(readmeLines(:), summaryXLSX, 'Sheet', 'README');
writecell(indexRows,      summaryXLSX, 'Sheet', 'Index');

writetable(excludedTable, summaryXLSX, 'Sheet', 'Excluded_Columns');
writetable(gapTable,      summaryXLSX, 'Sheet', 'GapFill_Report');
writetable(anchorTable,   summaryXLSX, 'Sheet', 'Anchor_Report');

writetable(T_full.Min360, summaryXLSX, 'Sheet', 'FullLadder_Min360');
writetable(T_full.Min60,  summaryXLSX, 'Sheet', 'FullLadder_Min60');
writetable(T_sel.Min360,  summaryXLSX, 'Sheet', 'Selective_Min360');
writetable(T_sel.Min60,   summaryXLSX, 'Sheet', 'Selective_Min60');
writetable(T_qc.Min360,   summaryXLSX, 'Sheet', 'QC_SelVsFull_Min360');
writetable(T_qc.Min60,    summaryXLSX, 'Sheet', 'QC_SelVsFull_Min60');
writetable(T_cross,       summaryXLSX, 'Sheet', 'QC_CrossMin_360vs60');

% LAST SHEET: Recommendation (explicitly written last)
writetable(recommendationTable, summaryXLSX, 'Sheet', 'Recommendation');

% Errors sheet
writetable(errorsTable, summaryXLSX, 'Sheet', 'Errors');

%% ------------------------- WRITE DETAIL WORKBOOK -------------------------

detailIndex = {
    'Sheet','Contents'
    'Index','Sheet map'
    'Sel_WindowLevel_Min360','Per-window fits used by Selective at Min360'
    'Sel_WindowLevel_Min60','Per-window fits used by Selective at Min60'
    'Suscept_Periods_360vs60','Wavelet global spectrum attenuation labels (FullLadder Min60 vs Min360)'
};
writecell(detailIndex, detailXLSX, 'Sheet', 'Index');
writetable(T_win.Min360, detailXLSX, 'Sheet', 'Sel_WindowLevel_Min360');
writetable(T_win.Min60,  detailXLSX, 'Sheet', 'Sel_WindowLevel_Min60');
writetable(T_sus,        detailXLSX, 'Sheet', 'Suscept_Periods_360vs60');

fprintf('\nDone.\n');
fprintf('Reports saved:\n  %s\n  %s\n', summaryXLSX, detailXLSX);

clc;
fprintf('Harmonic subtraction complete.\n');

%% -------------------------------------------------------------------------
% Local functions
% -------------------------------------------------------------------------
function s = tf(x)
    if x, s = 'OK'; else, s = 'NOT OK'; end
end

function v = light_value_repr(lightVec)
    % Return a compact representation of the Light duration condition (h).
    % If constant numeric, returns that scalar. Otherwise returns NaN.
    try
        if isnumeric(lightVec) || islogical(lightVec)
            u = unique(lightVec(~isnan(lightVec)));
            if isscalar(u)
                v = u;
            else
                v = NaN;
            end
        else
            v = NaN;
        end
    catch
        v = NaN;
    end
end

function [timeMinutes, TsMinutes] = infer_time_minutes(timeCol, timeName)
    if isdatetime(timeCol)
        t = timeCol(:);
        dt = minutes(diff(t));
        TsMinutes = median(dt, 'omitnan');
        timeMinutes = minutes(t - t(1));
        return;
    end

    if isnumeric(timeCol)
        t = timeCol(:);
        nm = lower(string(timeName));

        if contains(nm, "min")
            timeMinutes = t;
        elseif contains(nm, "hr") || contains(nm, "hour")
            timeMinutes = t * 60;
        elseif contains(nm, "day")
            timeMinutes = t * 24 * 60;
        else
            dtRaw = diff(t);
            medStep = median(dtRaw, 'omitnan');
            if max(t, [], 'omitnan') > 24 && medStep > 0 && medStep < 1
                timeMinutes = t * 60;
            else
                timeMinutes = t;
            end
        end

        dt = diff(timeMinutes);
        TsMinutes = median(dt, 'omitnan');
        if ~isfinite(TsMinutes) || TsMinutes <= 0
            error('Sampling interval could not be inferred from the Time column.');
        end
        return;
    end

    error('Unsupported Time column type. Time must be numeric or datetime.');
end

function [P0_h, maxDeltaR2, peakZ, Amp, ampSNR, deltaCurve] = regression_anchor_scan(t_hr, y, periods_h)
    t_hr = t_hr(:);
    y    = y(:);

    [deltaCurve, ampCurve, y0] = anchor_delta_curve(t_hr, y, periods_h);

    if all(~isfinite(deltaCurve))
        P0_h = NaN; maxDeltaR2 = NaN; peakZ = NaN; Amp = NaN; ampSNR = NaN;
        return;
    end

    [maxDeltaR2, iBest] = max(deltaCurve, [], 'omitnan');
    if ~isfinite(maxDeltaR2) || isempty(iBest) || iBest < 1 || iBest > numel(periods_h)
        P0_h = NaN; Amp = NaN;
    else
        P0_h = periods_h(iBest);
        Amp  = ampCurve(iBest);
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
    if ~isfinite(robustStd) || robustStd <= 0 || ~isfinite(Amp)
        ampSNR = 0;
    else
        ampSNR = Amp / robustStd;
    end
end

function bestDelta = regression_anchor_maxDeltaR2(t_hr, y, periods_h)
    t_hr = t_hr(:);
    y    = y(:);

    [deltaCurve, ~, ~] = anchor_delta_curve(t_hr, y, periods_h);
    if all(~isfinite(deltaCurve))
        bestDelta = NaN;
    else
        bestDelta = max(deltaCurve, [], 'omitnan');
        if ~isfinite(bestDelta), bestDelta = NaN; end
    end
end

function [deltaCurve, ampCurve, y0] = anchor_delta_curve(t_hr, y, periods_h)
    t_hr = t_hr(:);
    y    = y(:);

    tz = zscore_local(t_hr);
    X0 = [ones(size(t_hr)), tz, tz.^2];
    b0 = X0 \ y;
    y0 = X0 * b0;

    ssTot  = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
    ssRes0 = sum((y - y0).^2, 'omitnan');
    R2_0   = 1 - (ssRes0 / (ssTot + eps));

    deltaCurve = NaN(size(periods_h));
    ampCurve   = NaN(size(periods_h));

    for i = 1:numel(periods_h)
        P = periods_h(i);
        w = 2*pi / P;
        X = [X0, cos(w*t_hr), sin(w*t_hr)];
        b = X \ y;
        yhat = X * b;

        ssRes = sum((y - yhat).^2, 'omitnan');
        R2 = 1 - (ssRes / (ssTot + eps));

        dR2 = R2 - R2_0;
        deltaCurve(i) = dR2;

        a = b(end-1); s = b(end);
        ampCurve(i) = hypot(a, s);
    end
end

function yS = block_shuffle_surrogate(y, blockLenSamples)
    y = y(:);
    n = numel(y);

    if n < 2*blockLenSamples
        yS = y(randperm(n));
        return;
    end

    nBlocks = floor(n / blockLenSamples);
    remN    = n - nBlocks*blockLenSamples;

    blocks = cell(nBlocks + (remN>0), 1);

    for b = 1:nBlocks
        i1 = (b-1)*blockLenSamples + 1;
        i2 = b*blockLenSamples;
        blocks{b} = y(i1:i2);
    end

    if remN > 0
        blocks{end} = y(nBlocks*blockLenSamples + 1:end);
    end

    perm = randperm(numel(blocks));
    yS = vertcat(blocks{perm});
end

function z = zscore_local(x)
    x = x(:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        z = x - mu;
    else
        z = (x - mu) / sd;
    end
end

function [X0, Xh] = build_drift_and_harmonics_design(t_hr, t_min, P0_min, kList)
    tz = zscore_local(t_hr(:));
    X0 = [ones(size(tz)), tz, tz.^2];

    t = t_min(:);
    n = numel(t);
    nK = numel(kList);

    Xh = zeros(n, 2*nK);
    for ii = 1:nK
        kk = kList(ii);
        w = 2*pi*kk / P0_min;
        col = 2*(ii-1) + 1;
        Xh(:, col)   = cos(w*t);
        Xh(:, col+1) = sin(w*t);
    end
end

function [A, phi, R2] = fit_sinusoid_amp_phase_R2_with_linear(t_min, y, period_min)
    t = t_min(:);
    y = y(:);

    if numel(t) < 10 || all(~isfinite(y))
        A = NaN; phi = NaN; R2 = NaN;
        return;
    end

    w = 2*pi / period_min;
    tZ = zscore_local(t);
    X = [cos(w*t), sin(w*t), ones(size(t)), tZ];
    b = X \ y;
    yhat = X*b;

    a = b(1); s = b(2);
    A = hypot(a, s);
    phi = atan2(-s, a);

    ssRes = sum((y - yhat).^2, 'omitnan');
    ssTot = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
    if ssTot <= 0
        R2 = NaN;
    else
        R2 = 1 - (ssRes / ssTot);
    end
end

function [starts, ends] = make_windows(tMin, tMax, winLenMin, stepMin)
    starts = (tMin:stepMin:(tMax - winLenMin))';
    ends = starts + winLenMin;
    if isempty(starts)
        starts = tMin;
        ends = tMax;
    end
end

function x = wrapToPi_local(x)
    x = mod(x + pi, 2*pi) - pi;
end

function segs = contiguous_segments(mask)
    mask = mask(:);
    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts ends];
end

function [filledGapsCount, longGapsCount, unfilledSamples] = quantify_gaps(missingOriginal, missingAfterFill)
    missOrig  = missingOriginal(:);
    missAfter = missingAfterFill(:);

    segsOrig = contiguous_segments(missOrig);
    filledGapsCount = 0;
    longGapsCount = 0;

    for i = 1:size(segsOrig, 1)
        i1 = segsOrig(i, 1);
        i2 = segsOrig(i, 2);

        if ~any(missAfter(i1:i2))
            filledGapsCount = filledGapsCount + 1;
        else
            longGapsCount = longGapsCount + 1;
        end
    end

    unfilledSamples = sum(missOrig & missAfter);
end

function L = longest_nan_run(missingMask)
    if ~any(missingMask)
        L = 0;
        return;
    end
    segs = contiguous_segments(missingMask(:));
    L = max(segs(:,2) - segs(:,1) + 1);
end

function c = column_to_cell(x)
    x = x(:);
    if isdatetime(x)
        c = num2cell(x);
    elseif isnumeric(x) || islogical(x)
        c = num2cell(x);
    elseif isstring(x) || ischar(x) || iscellstr(x)
        c = cellstr(x);
    elseif iscell(x)
        c = x;
    else
        c = cellstr(string(x));
    end
end

function s = sanitise_filename(strIn)
    s = regexprep(string(strIn), '[^\w\-]', '_');
    s = char(s);
end

function [rFull, remFull] = map_back_to_full(x0, validIdx, residual, removed)
    rFull   = NaN(size(x0));
    remFull = NaN(size(x0));
    rFull(validIdx)   = residual;
    remFull(validIdx) = removed;

    rFull(isnan(x0))   = NaN;
    remFull(isnan(x0)) = NaN;
end

function v = variance_explained(y, residual)
    vy = var(y, 0, 'omitnan');
    if ~isfinite(vy) || vy <= 0
        v = NaN;
        return;
    end
    v = 1 - var(residual, 0, 'omitnan') / (vy + eps);
end

function write_output_excel(outFolder, fileName, timeName, timeOut, lightName, lightOut, dataIdx, varNames, dataMatrix, prefix)
    N = size(dataMatrix, 1);

    headers = cell(1, 2 + numel(dataIdx));
    headers{1} = timeName;
    for c = 1:numel(dataIdx)
        headers{1 + c} = [prefix '_' varNames{dataIdx(c)}];
    end
    headers{end} = lightName;

    out = cell(N+1, numel(headers));
    out(1,:) = headers;
    out(2:end,1) = column_to_cell(timeOut);
    out(2:end,2:(1+numel(dataIdx))) = num2cell(dataMatrix);
    out(2:end,end) = column_to_cell(lightOut);

    outFile = fullfile(outFolder, [prefix '_' fileName]);
    writecell(out, outFile);
    fprintf('%s saved: %s\n', prefix, outFile);
end