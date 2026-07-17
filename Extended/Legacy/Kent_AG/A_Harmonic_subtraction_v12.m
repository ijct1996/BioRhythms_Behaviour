% -------------------------------------------------------------------------
% Script: Harmonic_Subtraction_v12
% -------------------------------------------------------------------------
% Developed in version: MATLAB R2025a
% -------------------------------------------------------------------------
% Developed by: Isaiah J. Ting
% -------------------------------------------------------------------------
% Last updated: 08.05.2026
% -------------------------------------------------------------------------
clearvars; close all; clc;

%% --------------------------- RUN MODE SELECTION --------------------------

runMode = questdlg( ...
    'Choose input mode:', ...
    'Harmonic subtraction: input mode', ...
    'Single file', 'Multiple files', 'Cancel', ...
    'Single file');

if isempty(runMode) || strcmpi(runMode, 'Cancel')
    fprintf('No mode selected. Exiting.\n');
    return;
end

%% ------------------------------ FILE SELECTION ---------------------------

[fileList, inputParent] = select_input_files(runMode); %#ok<NASGU>
if isempty(fileList)
    fprintf('No input file(s) selected. Exiting.\n');
    return;
end

nFiles = numel(fileList);
fprintf('Selected %d file(s).\n', nFiles);
for iF = 1:nFiles
    fprintf('  %d) %s\n', iF, fileList{iF});
end

%% ------------------------------ OUTPUT ROOT -----------------------------

outputRoot = uigetdir(pwd, 'Select output root folder');
if isequal(outputRoot, 0)
    fprintf('No output folder selected. Exiting.\n');
    return;
end
fprintf('Output root: %s\n', outputRoot);

warn_if_path_long(outputRoot, 'Output root');

%% ------------------------- GLOBAL THRESHOLDS (ONCE) ---------------------

thresholds = ask_selective_thresholds_dialog();
fprintf('Thresholds: ratioTol=%.3f, Rcrit=%.3f, rhoCrit=%.3f, ampFracCritShort=%.3f\n', ...
    thresholds.ratioTol, thresholds.Rcrit, thresholds.rhoCrit, thresholds.ampFracCritShort);

%% -------------------- PREFLIGHT COLUMN MAPPING FOR ALL FILES ------------

fprintf('\nPreflight mapping (all files)...\n');

jobs = repmat(struct(), nFiles, 1);

for iF = 1:nFiles
    inputFile = fileList{iF};
    [~, fileStem, ext] = fileparts(inputFile);
    fileBaseWithExt = [fileStem ext];

    fprintf('\nPreflight %d/%d: %s\n', iF, nFiles, fileBaseWithExt);

    try
        [dataTablePF, readInfoPF] = read_input_table_preserve_robust(inputFile);
        [dataTablePF, droppedInfoPF] = drop_empty_columns_robust(dataTablePF);

        if isempty(dataTablePF) || width(dataTablePF) < 3
            error('Insufficient columns after import cleanup (need Time, >=1 Data, Light duration (h)).');
        end

        mapPF = preflight_column_mapping_dialog(dataTablePF, iF, nFiles);

        jobs(iF).inputFile = inputFile;
        jobs(iF).fileBaseWithExt = fileBaseWithExt;
        jobs(iF).fileStem = fileStem;
        jobs(iF).preflightMap = mapPF;
        jobs(iF).preflightReadInfo = readInfoPF;
        jobs(iF).preflightDroppedInfo = droppedInfoPF;

        if strcmpi(runMode, 'Multiple files')
            jobs(iF).outputFolder = unique_subfolder_for_file(outputRoot, iF, fileStem);
        else
            jobs(iF).outputFolder = outputRoot;
        end

        jobs(iF).thresholds = thresholds;

    catch ME
        error('Preflight failed for "%s": %s', inputFile, ME.message);
    end
end

fprintf('\nPreflight complete. Running analysis (no further mapping dialogs)...\n');

%% -------------------------- RUN ALL FILES -------------------------------

runErrors = cell(0,2);

for iF = 1:nFiles
    job = jobs(iF);

    fprintf('\n============================================================\n');
    fprintf('Running file %d/%d: %s\n', iF, nFiles, job.fileBaseWithExt);
    fprintf('Output: %s\n', job.outputFolder);
    fprintf('============================================================\n');

    try
        if ~exist(job.outputFolder, 'dir')
            mkdir(job.outputFolder);
        end

        warn_if_path_long(job.outputFolder, sprintf('Output folder (file %d)', iF));

        run_harmonic_subtraction_core(job);

    catch ME
        fprintf('FAILED: %s\n', ME.message);
        runErrors(end+1, :) = {job.fileBaseWithExt, ME.message}; %#ok<SAGROW>
    end
end

%% ----------------------------- FINAL SUMMARY ----------------------------

fprintf('\nBatch run complete.\n');
if isempty(runErrors)
    fprintf('All files processed successfully.\n');
else
    fprintf('Files failed (%d):\n', size(runErrors,1));
    for i = 1:size(runErrors,1)
        fprintf('  - %s | %s\n', runErrors{i,1}, runErrors{i,2});
    end
end

fprintf('Harmonic subtraction complete.\n');

%% ------------------------------------------------------------------------
% Core runner (single file, non-interactive)
% ------------------------------------------------------------------------
function run_harmonic_subtraction_core(job)

PIPELINE_VERSION = 'v12.0';
PIPELINE_NAME    = 'Harmonic_Subtraction_v12_validation_manifest';

inputFile    = job.inputFile;
outputFolder = job.outputFolder;
pfMap        = job.preflightMap;
thr          = job.thresholds;

[~, fileStem, ext] = fileparts(inputFile);
inputBaseWithExt = [fileStem ext];

%% NAMING CONVENTIONS
% Compact names reduce Windows path-length risk. Edit here to change naming globally.
NC = hs_naming();
SH = hs_sheets();

%% ------------------------------ SUBFOLDERS ------------------------------

reportsFolder = fullfile(outputFolder, NC.REPORTS_DIR);

tsRoot         = fullfile(outputFolder, NC.TS_DIR);
tsResidualFolder = fullfile(tsRoot, NC.TS_RES_DIR);
tsRemovedFolder  = fullfile(tsRoot, NC.TS_REM_DIR);

figAnchorFolder = fullfile(outputFolder, NC.FIG_ANCHOR_DIR);

scaloRoot       = fullfile(outputFolder, NC.SCALO_DIR);
scaloResFolder  = fullfile(scaloRoot, NC.SCALO_RES_DIR);
scaloRemFolder  = fullfile(scaloRoot, NC.SCALO_REM_DIR);

if ~exist(reportsFolder, 'dir');      mkdir(reportsFolder);      end
if ~exist(tsRoot, 'dir');             mkdir(tsRoot);             end
if ~exist(tsResidualFolder, 'dir');   mkdir(tsResidualFolder);   end
if ~exist(tsRemovedFolder,  'dir');   mkdir(tsRemovedFolder);    end
if ~exist(figAnchorFolder, 'dir');    mkdir(figAnchorFolder);    end

summaryXLSX = fullfile(reportsFolder, NC.SUMMARY_XLSX_NAME);
detailXLSX  = fullfile(reportsFolder, NC.DETAIL_XLSX_NAME);

enforce_output_path_length(summaryXLSX, 'Summary workbook');
enforce_output_path_length(detailXLSX,  'Detail workbook');

if exist(summaryXLSX, 'file'); delete(summaryXLSX); end
if exist(detailXLSX,  'file'); delete(detailXLSX);  end

%% ------------------------- READ DATA (ROBUST) ---------------------------

try
    [dataTable, ~] = read_input_table_preserve_robust(inputFile);
    [dataTable, droppedInfoNow] = drop_empty_columns_robust(dataTable);
catch ME
    error('Error reading file: %s', ME.message);
end

if isempty(dataTable) || width(dataTable) < 3
    error('Insufficient columns (need Time, >=1 Data, Light duration (h)).');
end

varNames = dataTable.Properties.VariableNames;

%% ---------------------- RESOLVE PREFLIGHT MAPPING -----------------------

[timeIdx, lightIdx, dataIdx, excludedNames, mappingNotes] = resolve_preflight_mapping(dataTable, pfMap);

if isempty(timeIdx) || isempty(lightIdx) || isempty(dataIdx)
    error('Resolved mapping is incomplete (Time/Light/Data).');
end

if ~isempty(mappingNotes)
    fprintf('Mapping notes:\n');
    for ii = 1:numel(mappingNotes)
        fprintf('  - %s\n', mappingNotes{ii});
    end
end
if ~isempty(droppedInfoNow.DroppedNames)
    fprintf('Auto-dropped empty columns: %s\n', strjoin(droppedInfoNow.DroppedNames, ', '));
end

timeName  = varNames{timeIdx};
lightName = varNames{lightIdx};

fprintf('\nResolved columns:\n');
fprintf('  Time: %s\n', timeName);
fprintf('  Light duration (h): %s\n', lightName);
fprintf('  Data: %d columns\n', numel(dataIdx));
if ~isempty(excludedNames)
    fprintf('  Excluded: %d columns\n', numel(excludedNames));
end

%% ------------------------- FIXED INTERNAL SPECS -------------------------

minPeriodsToRun = [360, 60];   % Min 360 then Min 60

missingFracThreshold = 0.05;
MaxGapMinutes        = 30;

anchorBand_h  = [22, 28];
anchorStep_h  = 0.01;

blockLenHours = 8;
Nsurr         = 300;
alphaAnchor   = 0.05;

minCyclesForAnchor = 3.5;
minDeltaR2         = 0.005;
useEdgeGuard       = true;
edgeMargin_h       = 0.10;

PeakZ_warn = 1.0;
SNR_warn   = 1.0;

peakSearchFrac   = 0.06;
ratioTol         = thr.ratioTol;
Rcrit            = thr.Rcrit;
rhoCrit          = thr.rhoCrit;
ampFracCritShort = thr.ampFracCritShort;

% v10: locked vs free evidence tolerance
dR2_freeMinusLocked_Tol = 0.02;

stepHours = 24;
minSamplesForAnalysis = 250;
minSamplesPerWindow   = 80;

SUBTRACT_FUNDAMENTAL_IN_SELECTIVE = true;

SAVE_ANCHOR_FIGS = true;

SAVE_SCALOGRAMS = true;
SAVE_SCALOGRAMS_SELECTIVE = true;

SCALO_MIN_PERIOD_MIN = 60;
SCALO_MAX_PERIOD_MIN = 1590;
SCALO_YTICKS_H = 0:4:26;

nRMS_low  = 0.10;
nRMS_high = 0.25;
dVar_low  = 0.05;
dVar_high = 0.15;

% v10: internal auto anchor routing
ANCHOR_MODE = 'auto'; % 'auto' | 'fixed' | 'timevarying' (no UI)
anchorNonstat_IQRFrac_Thresh = 0.015;
anchorNonstat_Slope_hPerDay_Thresh = 0.10;
minWindowsForAnchorNonstat = 3;

%% --------------------- INFER TIME BASE AND SAMPLING ---------------------

timeCol = dataTable{:, timeIdx};
try
    [timeMinutesAll, TsMinutes] = infer_time_minutes(timeCol, timeName);
catch ME
    error('Error inferring time base: %s', ME.message);
end

if ~isfinite(TsMinutes) || TsMinutes <= 0
    error('Invalid sampling interval. Exiting.');
end

N = height(dataTable);
maxGapSamples = max(0, round(MaxGapMinutes / TsMinutes));

durationMinAll = max(timeMinutesAll, [], 'omitnan') - min(timeMinutesAll, [], 'omitnan');
durationHoursAll = durationMinAll / 60;

baselineWinHours = min(max(72, 0.5*durationHoursAll), 168);
baselineMinPoints = 20;
baselineWinSamples = max(baselineMinPoints, round((baselineWinHours*60) / TsMinutes));

blockLenSamples = max(4, round((blockLenHours*60) / TsMinutes));

fprintf('\nSampling: %.6g min/sample\n', TsMinutes);
fprintf('Duration: %.2f h\n', durationHoursAll);
fprintf('Baseline removal: movmedian %.1f h (%d samples)\n', baselineWinHours, baselineWinSamples);
fprintf('Anchor band: %.1f to %.1f h (step %.3f h)\n', anchorBand_h(1), anchorBand_h(2), anchorStep_h);
fprintf('Block-shuffle: %d surrogates | block %.1f h (%d samples)\n', Nsurr, blockLenHours, blockLenSamples);
fprintf('Runs: MinPeriods = [%s] min\n', num2str(minPeriodsToRun));

%% ------------------------- LIGHT CONDITION CHANGES ----------------------

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

time_min = timeMinutesAll(:);
time_day = time_min / (60*24);

%% ---------------------------- REPORT TABLES -----------------------------

nCols = numel(dataIdx);

gapReportHeader = {'LightDur_h','Column','MissingFraction','Interpolated','FilledSamples','FilledGapsCount','LongGapsCount','UnfilledSamples','LongestGapSamples'};
gapReport = cell(nCols, numel(gapReportHeader));
gap_i = 0;

anchorHeader = {'LightDur_h','Column','AnchorOK','Period_hours','MaxDeltaR2','PeakZ','AnchorAmp','AmpSNR','Duration_h','CyclesAtPeriod','pBlock', ...
                'AnchorModeUsed','P0w_median_h','P0w_IQR_h','P0w_IQRFrac','P0w_slope_hPerDay','AnchorNonstatIndex','Notes'};
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

selHeader  = {'LightDur_h','MinPeriod_mins','Column','AnchorOK','Period_hours','AnchorModeUsed','RemovedFundamental','K_removed','Harmonics_k','VarExplained','NumWindowsUsed', ...
              'kLikely','KLikely_max','Notes'};
fullHeader = {'LightDur_h','MinPeriod_mins','Column','AnchorOK','Period_hours','AnchorModeUsed','K_removed','VarExplained','Notes'};

selReport = struct();
fullReport = struct();
selWinReport = struct();

selWinHeader = {'LightDur_h','MinPeriod_mins','Column','WindowIndex','WinStart_h','WinEnd_h','Period_hours','k', ...
                'Amp_locked','Phase_locked_rad','FitR2_locked','FitR2_freeBest','dR2_freeMinusLocked','PkBest_h','Ratio','Note'};

qcSelFullHeader = {'LightDur_h','MinPeriod_mins','Column','AnchorOK','Period_hours','AnchorModeUsed', ...
    'VarExpl_Sel','VarExpl_Full','DeltaVarExpl_FullMinusSel','RMSdiff','nRMS','Risk','Notes'};
qcSelFull = struct();

crossMinHeader = {'LightDur_h','Column','AnchorOK','Period_hours','AnchorModeUsed', ...
    'K_Min360','K_Min60','VarExplFull_Min360','VarExplFull_Min60','DeltaVarExpl_60minus360','RMS_60vs360','nRMS_60vs360', ...
    'KLikelyMax_Min360','KLikelyMax_Min60','Recommendation','Notes'};
crossMinQC = cell(nCols, numel(crossMinHeader));
cross_i = 0;

susceptHeader = {'LightDur_h','Column','Period_min','GWS_Min360','GWS_Min60','AttenuationRatio','Label'};
maxSusRows = max(1, nCols * 250);
susceptTable = cell(maxSusRows, numel(susceptHeader));
suscept_i = 0;

errors = cell(nCols, 2);
err_i  = 0;

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

%% ----------------------- OUTPUT ARRAYS (FULL LENGTH) --------------------

resFull = struct(); remFull = struct();
resSel  = struct(); remSel  = struct();

for mp = minPeriodsToRun
    key = sprintf('Min%d', mp);
    resFull.(key) = NaN(N, numel(dataIdx));
    remFull.(key) = NaN(N, numel(dataIdx));
    resSel.(key)  = NaN(N, numel(dataIdx));
    remSel.(key)  = NaN(N, numel(dataIdx));
end

%% ----------------------- PROCESS EACH DATA COLUMN -----------------------

for c = 1:numel(dataIdx)
    colName = varNames{dataIdx(c)};
    fprintf('\nProcessing %d/%d: %s\n', c, numel(dataIdx), colName);

    try
        x0 = dataTable{:, dataIdx(c)};
        if ~isnumeric(x0)
            x0 = str2double(string(x0));
        end
        x0 = x0(:);

        LightDur_h = light_value_repr(lightVec);

        % Missingness and gap fill
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

        % Baseline detrend
        base = movmedian(xFill, baselineWinSamples, 'omitnan', 'Endpoints', 'shrink');
        yDetrFull = xFill - base;

        validIdx = find(~isnan(yDetrFull) & isfinite(timeMinutesAll));
        if numel(validIdx) < minSamplesForAnalysis
            fprintf('  Too few valid samples; skipping subtraction for this column.\n');

            for mp = minPeriodsToRun
                key = sprintf('Min%d', mp);

                resFull.(key)(:, c) = x0;  remFull.(key)(:, c) = NaN(size(x0));
                resSel.(key)(:, c)  = x0;  remSel.(key)(:, c)  = NaN(size(x0));

                full_i.(key) = full_i.(key) + 1;
                fullReport.(key)(full_i.(key), :) = {LightDur_h, mp, colName, false, NaN, 'NA', 0, 0, 'Too few valid samples'};

                sel_i.(key) = sel_i.(key) + 1;
                selReport.(key)(sel_i.(key), :)  = {LightDur_h, mp, colName, false, NaN, 'NA', SUBTRACT_FUNDAMENTAL_IN_SELECTIVE, 0, '', 0, 0, '', NaN, 'Too few valid samples'};

                qc_i.(key) = qc_i.(key) + 1;
                qcSelFull.(key)(qc_i.(key), :)  = {LightDur_h, mp, colName, false, NaN, 'NA', NaN, NaN, NaN, NaN, NaN, 'NA', 'Too few valid samples'};
            end

            anchor_i = anchor_i + 1;
            anchorReport(anchor_i, :) = {LightDur_h, colName, false, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                                         'NA', NaN, NaN, NaN, NaN, NaN, 'Too few valid samples'};
            continue;
        end

        t_min = timeMinutesAll(validIdx);
        t_hr  = t_min / 60;

        durCol_h = (max(t_min, [], 'omitnan') - min(t_min, [], 'omitnan')) / 60;

        y = yDetrFull(validIdx);
        y = y - median(y, 'omitnan');

        % Anchor scan
        periods_h = (anchorBand_h(1):anchorStep_h:anchorBand_h(2))';
        [P0_h, maxDeltaR2, peakZ, A0, ampSNR, deltaCurve] = regression_anchor_scan(t_hr, y, periods_h);
        cyclesAtP0 = durCol_h / max(P0_h, eps);

        % Block-shuffle p-value
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

        % Anchor gating
        reasons = {};
        warnings = {};

        finiteOK = isfinite(P0_h) && isfinite(maxDeltaR2) && isfinite(pBlock) && isfinite(cyclesAtP0);
        if ~finiteOK, reasons = [reasons, {'Non-finite anchor statistics'}]; end

        sigOK = (pBlock < alphaAnchor);
        if ~sigOK, reasons = [reasons, {sprintf('pBlock too large (%.3g >= %.3g)', pBlock, alphaAnchor)}]; end

        cyclesOK = (cyclesAtP0 >= minCyclesForAnchor);
        if ~cyclesOK, reasons = [reasons, {sprintf('Too few cycles (%.2f < %.1f)', cyclesAtP0, minCyclesForAnchor)}]; end

        deltaOK = (maxDeltaR2 >= minDeltaR2);
        if ~deltaOK, reasons = [reasons, {sprintf('DeltaR2 too small (%.4g < %.4g)', maxDeltaR2, minDeltaR2)}]; end

        edgeOK = true;
        if useEdgeGuard
            edgeOK = (P0_h > (anchorBand_h(1) + edgeMargin_h)) && (P0_h < (anchorBand_h(2) - edgeMargin_h));
            if ~edgeOK, reasons = [reasons, {sprintf('Period near band edge (P=%.3f h)', P0_h)}]; end
        end

        if isfinite(peakZ) && peakZ < PeakZ_warn
            warnings = [warnings, {sprintf('QC: broad peak (PeakZ=%.2f < %.2f)', peakZ, PeakZ_warn)}]; %#ok<AGROW>
        end
        if isfinite(ampSNR) && ampSNR < SNR_warn
            warnings = [warnings, {sprintf('QC: low SNR (AmpSNR=%.2f < %.2f)', ampSNR, SNR_warn)}];
        end

        anchorOK = finiteOK && sigOK && cyclesOK && deltaOK && edgeOK;

        if isempty(reasons)
            noteAnchor = 'Anchor accepted';
        else
            noteAnchor = strjoin(reasons, '; ');
        end
        if ~isempty(warnings)
            noteAnchor = [noteAnchor ' | ' strjoin(warnings, '; ')];
        end

        fprintf('  Anchor: %s | P=%.3f h | dR2=%.4g | Z=%.2f | Amp=%.4g | SNR=%.2f | cycles=%.2f | p=%.3g\n', ...
            tf(anchorOK), P0_h, maxDeltaR2, peakZ, A0, ampSNR, cyclesAtP0, pBlock);

        % Anchor nonstationarity diagnostics and auto routing
        anchorModeUsed = 'fixed';
        P0w_median = NaN; P0w_IQR = NaN; P0w_IQRFrac = NaN; P0w_slope = NaN; anchorNonstatIndex = NaN;

        winHours_diag = min(max(72, round(durCol_h/2)), 168);
        [winStarts_min_diag, winEnds_min_diag] = make_windows(min(t_min), max(t_min), winHours_diag*60, stepHours*60);
        [P0w_vec, tCenter_h_vec, nWinUsed] = windowed_anchor_periods(t_min, y, winStarts_min_diag, winEnds_min_diag, periods_h, minSamplesPerWindow);

        if nWinUsed >= minWindowsForAnchorNonstat
            P0w_median = median(P0w_vec, 'omitnan');
            P0w_IQR    = iqr(P0w_vec);
            P0w_IQRFrac = P0w_IQR / (P0w_median + eps);

            try
                tt = tCenter_h_vec(:);
                pp = P0w_vec(:);
                ok = isfinite(tt) & isfinite(pp);
                if sum(ok) >= 3
                    b = polyfit(tt(ok), pp(ok), 1); % slope per hour
                    P0w_slope = b(1) * 24;         % slope per day
                end
            catch
                P0w_slope = NaN;
            end

            anchorNonstatIndex = max([P0w_IQRFrac, abs(P0w_slope)/(P0w_median+eps)], [], 'omitnan');
        end

        if strcmpi(ANCHOR_MODE, 'fixed')
            anchorModeUsed = 'fixed';
        elseif strcmpi(ANCHOR_MODE, 'timevarying')
            anchorModeUsed = 'timevarying';
        else
            useTV = false;
            if isfinite(P0w_IQRFrac) && (P0w_IQRFrac >= anchorNonstat_IQRFrac_Thresh)
                useTV = true;
            end
            if isfinite(P0w_slope) && (abs(P0w_slope) >= anchorNonstat_Slope_hPerDay_Thresh)
                useTV = true;
            end
            if nWinUsed < minWindowsForAnchorNonstat
                useTV = false;
            end
            anchorModeUsed = ternary(useTV, 'timevarying', 'fixed');
        end

        anchor_i = anchor_i + 1;
        anchorReport(anchor_i, :) = {LightDur_h, colName, anchorOK, P0_h, maxDeltaR2, peakZ, A0, ampSNR, durCol_h, cyclesAtP0, pBlock, ...
                                     anchorModeUsed, P0w_median, P0w_IQR, P0w_IQRFrac, P0w_slope, anchorNonstatIndex, noteAnchor};

        % Anchor figure
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
                title(sprintf('%s | P %.2f h', colName, P0_h), 'FontName','Times New Roman', 'Interpreter','none');

                nexttile;
                histogram(surrMax, 25);
                hold on; xline(maxDeltaR2, 'LineWidth', 1.2); hold off;
                set(gca, 'Box','off', 'FontName','Times New Roman', 'TickDir','out', ...
                    'XGrid','off','YGrid','off');
                xlabel('Surrogate max DeltaR^2 (22 to 28 h)', 'FontWeight','bold', 'FontName','Times New Roman');
                ylabel('Count', 'FontWeight','bold', 'FontName','Times New Roman');
                title(sprintf('Block-shuffle null | p=%.3g', pBlock), 'FontName','Times New Roman');

                outFig = fullfile(figAnchorFolder, sprintf('%s_%s.jpg', NC.ANCH_PREFIX, sanitise_filename(colName)));
                print(fig, outFig, '-djpeg', '-r600');
                close(fig);
            catch
            end
        end

        % If anchor rejected: pass-through
        if ~anchorOK
            for mp = minPeriodsToRun
                key = sprintf('Min%d', mp);

                resFull.(key)(:, c) = x0;  remFull.(key)(:, c) = NaN(size(x0));
                resSel.(key)(:, c)  = x0;  remSel.(key)(:, c)  = NaN(size(x0));

                full_i.(key) = full_i.(key) + 1;
                fullReport.(key)(full_i.(key), :) = {LightDur_h, mp, colName, false, P0_h, anchorModeUsed, 0, 0, 'Anchor rejected (no subtraction)'};

                sel_i.(key) = sel_i.(key) + 1;
                selReport.(key)(sel_i.(key), :)  = {LightDur_h, mp, colName, false, P0_h, anchorModeUsed, SUBTRACT_FUNDAMENTAL_IN_SELECTIVE, 0, '', 0, 0, '', NaN, 'Anchor rejected (no subtraction)'};

                qc_i.(key) = qc_i.(key) + 1;
                qcSelFull.(key)(qc_i.(key), :)  = {LightDur_h, mp, colName, false, P0_h, anchorModeUsed, NaN, NaN, NaN, NaN, NaN, 'NA', 'Anchor rejected'};
            end
            continue;
        end

        % Anchor phase basis
        if strcmpi(anchorModeUsed, 'fixed')
            phi0 = 2*pi * (t_hr / P0_h);
        else
            try
                [phi0, ~] = build_timevarying_anchor_phase(t_hr, t_min, y, winStarts_min_diag, winEnds_min_diag, periods_h, minSamplesPerWindow, P0_h);
            catch
                phi0 = 2*pi * (t_hr / P0_h);
                anchorReport{anchor_i, end} = [anchorReport{anchor_i, end} ' | TV anchor failed: fallback to fixed'];
                anchorModeUsed = 'fixed';
            end
        end

        % Windowing for coupling
        winHours = min(max(72, round(durCol_h/2)), 168);

        tMinAll = min(t_min); tMaxAll = max(t_min);
        [winStarts_min, winEnds_min] = make_windows(tMinAll, tMaxAll, winHours*60, stepHours*60);
        nW = numel(winStarts_min);

        fprintf('  Windows: %g h | step %g h | nW=%d | AnchorMode=%s\n', winHours, stepHours, nW, anchorModeUsed);

        % Fundamental per-window fit (locked to phi0)
        A1_w   = NaN(nW,1);
        Phi1_w = NaN(nW,1);
        P0_w_forTarget = NaN(nW,1);

        for w = 1:nW
            wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
            if sum(wMask) < minSamplesPerWindow
                continue;
            end

            tmw_hr = t_hr(wMask);
            ymw    = y(wMask);
            phiw   = phi0(wMask);

            [A1, phi1, R2w] = fit_locked_phase_amp_phase_R2_with_linear(phiw, tmw_hr, ymw);

            A1_w(w)   = A1;
            Phi1_w(w) = phi1;

            try
                P0w = window_anchor_period_single(tmw_hr, ymw, periods_h);
            catch
                P0w = P0_h;
            end
            P0_w_forTarget(w) = P0w;

            for mp = minPeriodsToRun
                key = sprintf('Min%d', mp);
                win_i.(key) = win_i.(key) + 1;
                selWinReport.(key)(win_i.(key), :) = {LightDur_h, mp, colName, w, winStarts_min(w)/60, winEnds_min(w)/60, ...
                    P0w, 1, A1, phi1, R2w, R2w, 0, P0w, 1, 'Fundamental (locked to anchor phase, +linear)'};
            end
        end

        okFund = isfinite(A1_w) & isfinite(Phi1_w) & isfinite(P0_w_forTarget);
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

        % Prepare scalogram filterbank
        if SAVE_SCALOGRAMS
            if ~exist(scaloRoot, 'dir');      mkdir(scaloRoot);      end
            if ~exist(scaloResFolder, 'dir'); mkdir(scaloResFolder); end
            if ~exist(scaloRemFolder, 'dir'); mkdir(scaloRemFolder); end

            scaloFB = cwtfilterbank('SignalLength', N, ...
                'SamplingPeriod', minutes(TsMinutes), ...
                'PeriodLimits', [minutes(SCALO_MIN_PERIOD_MIN), minutes(SCALO_MAX_PERIOD_MIN)], ...
                'Wavelet', 'amor');
        end

        tmpFull = struct();
        tmpVarExplFull = struct();
        tmpK = struct();
        tmpKLikelyMax = struct();

        for mp = minPeriodsToRun
            key = sprintf('Min%d', mp);
            minPeriodInterestMin = mp;

            P0_min = P0_h * 60;
            K = floor(P0_min / minPeriodInterestMin);
            K = max(2, min(K, 24));
            kList = 2:K;

            Ak_lock_w     = NaN(nW, numel(kList));
            Phik_lock_w   = NaN(nW, numel(kList));
            R2_lock_w     = NaN(nW, numel(kList));

            R2_free_w     = NaN(nW, numel(kList));
            ratio_w       = NaN(nW, numel(kList));
            Pbest_h       = NaN(nW, numel(kList));
            dR2_freeMinusLocked_w = NaN(nW, numel(kList));

            for w = 1:nW
                wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
                if sum(wMask) < minSamplesPerWindow || ~isfinite(A1_w(w)) || ~isfinite(Phi1_w(w)) || ~isfinite(P0_w_forTarget(w))
                    continue;
                end

                tmw_min = t_min(wMask);
                tmw_hr  = t_hr(wMask);
                ymw     = y(wMask);
                phiw    = phi0(wMask);

                for j = 1:numel(kList)
                    k = kList(j);

                    targetP_min = (P0_w_forTarget(w) * 60) / k;
                    scanP_min = linspace(targetP_min*(1-peakSearchFrac), targetP_min*(1+peakSearchFrac), 19);

                    [AkL, phiL, R2L] = fit_locked_phase_amp_phase_R2_with_linear(k*phiw, tmw_hr, ymw);

                    bestR2  = -Inf;
                    bestP   = NaN;
                    for sp = 1:numel(scanP_min)
                        [~, ~, R2kk] = fit_sinusoid_amp_phase_R2_with_linear(tmw_min, ymw, scanP_min(sp));
                        if isfinite(R2kk) && R2kk > bestR2
                            bestR2 = R2kk;
                            bestP  = scanP_min(sp);
                        end
                    end

                    if ~isfinite(R2L)
                        continue;
                    end
                    if ~isfinite(bestR2)
                        bestR2 = NaN;
                    end

                    Ak_lock_w(w, j)   = AkL;
                    Phik_lock_w(w, j) = phiL;
                    R2_lock_w(w, j)   = R2L;

                    R2_free_w(w, j)   = bestR2;

                    if isfinite(bestP)
                        ratio_w(w, j) = bestP / targetP_min;
                        Pbest_h(w, j) = bestP / 60;
                    else
                        ratio_w(w, j) = NaN;
                        Pbest_h(w, j) = NaN;
                    end

                    dR2_freeMinusLocked_w(w, j) = bestR2 - R2L;

                    win_i.(key) = win_i.(key) + 1;
                    selWinReport.(key)(win_i.(key), :) = {LightDur_h, mp, colName, w, winStarts_min(w)/60, winEnds_min(w)/60, ...
                        P0_w_forTarget(w), k, AkL, phiL, R2L, bestR2, dR2_freeMinusLocked_w(w, j), Pbest_h(w, j), ratio_w(w, j), ...
                        'Harmonic candidate (locked phase + free scan evidence, +linear)'};
                end
            end

            harmonicLikely = false(1, numel(kList));

            for j = 1:numel(kList)
                k = kList(j);

                okW = isfinite(A1_w) & isfinite(Ak_lock_w(:, j)) & isfinite(Phi1_w) & isfinite(Phik_lock_w(:, j)) & isfinite(ratio_w(:, j));
                nOk = sum(okW);

                if nOk < minWindowsForDecision
                    continue;
                end

                ratios = ratio_w(okW, j);
                ratioMad = mad(ratios - 1, 1);
                ratioPass = ratioMad <= ratioTol;

                del = wrapToPi_local(Phik_lock_w(okW, j) - k * Phi1_w(okW));
                Rplv = abs(mean(exp(1i * del)));
                Rpass = isfinite(Rplv) && (Rplv >= Rcrit);

                dR2v = dR2_freeMinusLocked_w(okW, j);
                dR2_med = median(dR2v, 'omitnan');
                lockedOK = isfinite(dR2_med) && (dR2_med <= dR2_freeMinusLocked_Tol);

                if strcmp(regime, 'long')
                    rho = corr(A1_w(okW), Ak_lock_w(okW, j), 'Type', 'Spearman', 'Rows', 'complete');
                    rhoPass = isfinite(rho) && (abs(rho) >= rhoCrit);
                    harmonicLikely(j) = ((ratioPass + rhoPass + Rpass) >= 2) && lockedOK;
                elseif strcmp(regime, 'short')
                    A1_med = median(A1_w(okW), 'omitnan');
                    Ak_med = median(Ak_lock_w(okW, j), 'omitnan');
                    ampPass = isfinite(A1_med) && A1_med > 0 && isfinite(Ak_med) && (Ak_med >= ampFracCritShort * A1_med);
                    harmonicLikely(j) = ratioPass && Rpass && ampPass && lockedOK;
                end
            end

            kLikely = kList(harmonicLikely);
            if isempty(kLikely)
                KLikely_max = NaN;
            else
                KLikely_max = max(kLikely);
            end
            tmpKLikelyMax.(key) = KLikely_max;

            % Selective subtraction
            if SUBTRACT_FUNDAMENTAL_IN_SELECTIVE
                kSubSel = unique([1, kLikely(:)'], 'stable');
            else
                kSubSel = kLikely(:)';
            end

            yFit = y;

            if isempty(kSubSel)
                removedSel = zeros(size(yFit));
                residualSel = yFit;
                varExplSel = 0;
                noteSel = sprintf('Anchor OK (P=%.3f h) but no harmonics classified (%s regime, nW=%d).', P0_h, regime, nFund);
            else
                [X0, Xh] = build_drift_and_harmonics_design_phi(t_hr, phi0, kSubSel);
                X = [X0, Xh];
                beta = X \ yFit;
                yhatDrift = X0 * beta(1:size(X0,2));
                yhatTotal = X  * beta;

                removedSel  = yhatTotal - yhatDrift;
                residualSel = yFit - removedSel;

                varExplSel = variance_explained(yFit, residualSel);
                noteSel = sprintf('Selective removed k=%s (%s regime, nW=%d).', mat2str(kSubSel), regime, nFund);
            end

            [rSelFull, remSelFull_] = map_back_to_full(x0, validIdx, residualSel, removedSel);
            resSel.(key)(:, c) = rSelFull;
            remSel.(key)(:, c) = remSelFull_;

            sel_i.(key) = sel_i.(key) + 1;
            selReport.(key)(sel_i.(key), :) = {LightDur_h, mp, colName, true, P0_h, anchorModeUsed, SUBTRACT_FUNDAMENTAL_IN_SELECTIVE, ...
                numel(kSubSel), mat2str(kSubSel), varExplSel, nFund, mat2str(kLikely), KLikely_max, noteSel};

            % Full ladder subtraction (v9 behaviour but with phase basis)
            kSubFull = 1:K;

            [X0F, XhF] = build_drift_and_harmonics_design_phi(t_hr, phi0, kSubFull);
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
            fullReport.(key)(full_i.(key), :) = {LightDur_h, mp, colName, true, P0_h, anchorModeUsed, K, varExplFullV, ...
                sprintf('Full ladder removed k=1..%d (minPeriod=%g min)', K, minPeriodInterestMin)};

            fprintf('  %s | Sel VarExpl=%.3f | FL VarExpl=%.3f (K=%d)\n', key, varExplSel, varExplFullV, K);

            % Selective vs Full ladder QC
            rmsDiff = sqrt(mean((residualSel(:) - residualFull(:)).^2, 'omitnan'));
            denom   = sqrt(mean((yFit(:)).^2, 'omitnan'));
            nRMS    = rmsDiff / (denom + eps);

            dVar = varExplFullV - varExplSel;

            risk = 'Low';
            if (nRMS >= nRMS_high) || (dVar >= dVar_high)
                risk = 'High';
            elseif (nRMS >= nRMS_low) || (dVar >= dVar_low)
                risk = 'Moderate';
            end

            qc_i.(key) = qc_i.(key) + 1;
            qcSelFull.(key)(qc_i.(key), :) = {LightDur_h, mp, colName, true, P0_h, anchorModeUsed, ...
                varExplSel, varExplFullV, dVar, rmsDiff, nRMS, risk, sprintf('Regime=%s, nW=%d, K=%d, KLikelyMax=%s', regime, nFund, K, num2str(KLikely_max))};

            % Scalograms
            if SAVE_SCALOGRAMS
                try
                    if ~exist(scaloRoot, 'dir');      mkdir(scaloRoot);      end
                    if ~exist(scaloResFolder, 'dir'); mkdir(scaloResFolder); end
                    if ~exist(scaloRemFolder, 'dir'); mkdir(scaloRemFolder); end

                    safeCol = sanitise_filename(colName);

                    % Labels for this min period
                    labelFL  = hs_label(NC.METHOD_FULL, mp, NC);
                    labelSEL = hs_label(NC.METHOD_SEL,  mp, NC);

                    % Full ladder residual/removed
                    sResF = rFullFull(:);    sResF(~isfinite(sResF)) = 0;
                    sRemF = remFullFull_(:); sRemF(~isfinite(sRemF)) = 0;

                    [wtResF, periodsResF] = cwt(sResF, 'FilterBank', scaloFB);
                    [wtRemF, periodsRemF] = cwt(sRemF, 'FilterBank', scaloFB);

                    pH_resF = hours(periodsResF);
                    pH_remF = hours(periodsRemF);

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
                    title(sprintf('Scalogram | Residual | FL %s | %s', key, colName), 'FontName','Times New Roman', 'Interpreter','none');
                    outResF = fullfile(scaloResFolder, hs_scalo_fig_name(labelFL, false, safeCol, NC));
                    print(figS1, outResF, '-djpeg', '-r600');
                    close(figS1);

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
                    title(sprintf('Scalogram | Removed | FL %s | %s', key, colName), 'FontName','Times New Roman', 'Interpreter','none');
                    outRemF = fullfile(scaloRemFolder, hs_scalo_fig_name(labelFL, true, safeCol, NC));
                    print(figS2, outRemF, '-djpeg', '-r600');
                    close(figS2);

                    if SAVE_SCALOGRAMS_SELECTIVE
                        sResS = rSelFull(:);      sResS(~isfinite(sResS)) = 0;
                        sRemS = remSelFull_(:);   sRemS(~isfinite(sRemS)) = 0;

                        [wtResS, periodsResS] = cwt(sResS, 'FilterBank', scaloFB);
                        [wtRemS, periodsRemS] = cwt(sRemS, 'FilterBank', scaloFB);

                        pH_resS = hours(periodsResS);
                        pH_remS = hours(periodsRemS);

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
                        title(sprintf('Scalogram | Residual | SEL %s | %s', key, colName), 'FontName','Times New Roman', 'Interpreter','none');
                        outResS = fullfile(scaloResFolder, hs_scalo_fig_name(labelSEL, false, safeCol, NC));
                        print(figSS1, outResS, '-djpeg', '-r600');
                        close(figSS1);

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
                        title(sprintf('Scalogram | Removed | SEL %s | %s', key, colName), 'FontName','Times New Roman', 'Interpreter','none');
                        outRemS = fullfile(scaloRemFolder, hs_scalo_fig_name(labelSEL, true, safeCol, NC));
                        print(figSS2, outRemS, '-djpeg', '-r600');
                        close(figSS2);
                    end
                catch
                end
            end

            tmpFull.(key) = residualFull(:);
            tmpVarExplFull.(key) = varExplFullV;
            tmpK.(key) = K;
        end

        % Cross-min QC
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
            crossMinQC(cross_i, :) = {LightDur_h, colName, true, P0_h, anchorModeUsed, ...
                tmpK.Min360, tmpK.Min60, tmpVarExplFull.Min360, tmpVarExplFull.Min60, dVar_cross, rmsCross, nRMS_cross, ...
                tmpKLikelyMax.Min360, tmpKLikelyMax.Min60, rec, ...
                sprintf('Thresholds: nRMS(%.2f/%.2f), dVar(%.2f/%.2f)', nRMS_low, nRMS_high, dVar_low, dVar_high)};
        end

        % Susceptibility periods (Full ladder only)
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
                        susceptTable = [susceptTable; cell(maxSusRows, numel(susceptHeader))]; %#ok<AGROW>
                    end
                    susceptTable(suscept_i, :) = {LightDur_h, colName, perMin(iP), gws360(iP), gws60(iP), ratio, label};
                end
            catch
            end
        end

    catch ME
        fprintf('  Error: %s\n', ME.message);
        err_i = err_i + 1;
        errors(err_i, :) = {inputBaseWithExt, sprintf('Column %s: %s', colName, ME.message)};
    end
end

%% ------------------------ TRIM PREALLOCATED LOGS -------------------------

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

%% -------------------- WRITE TIME-SERIES WORKBOOKS -----------------------

try
    timeOut  = dataTable{:, timeIdx};
    lightOut = dataTable{:, lightIdx};

    for mp = minPeriodsToRun
        key = sprintf('Min%d', mp); %#ok<NASGU>

        lblFL  = hs_label(NC.METHOD_FULL, mp, NC);
        lblSEL = hs_label(NC.METHOD_SEL,  mp, NC);

        write_output_excel(tsResidualFolder, fileStem, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            resFull.(sprintf('Min%d', mp)), lblFL);

        write_output_excel(tsRemovedFolder, fileStem, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            remFull.(sprintf('Min%d', mp)), lblFL);

        write_output_excel(tsResidualFolder, fileStem, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            resSel.(sprintf('Min%d', mp)), lblSEL);

        write_output_excel(tsRemovedFolder, fileStem, timeName, timeOut, lightName, lightOut, dataIdx, varNames, ...
            remSel.(sprintf('Min%d', mp)), lblSEL);
    end
catch ME
    fprintf('Output write failed: %s\n', ME.message);
end

%% ------------------------ WRITE REPORT WORKBOOKS ------------------------

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

T_cross = cell2table(crossMinQC,   'VariableNames', crossMinHeader);
T_sus   = cell2table(susceptTable, 'VariableNames', susceptHeader);

%% ------------------------ BUILD RECOMMENDATION SHEET ---------------------

recHeader = {'LightDur_h','Column','AnchorOK','Period_hours','AnchorModeUsed','RecommendedResidual','Workbook','SignalColumnHeader','Reason'};
recRows = cell(0, numel(recHeader));

for c = 1:numel(dataIdx)
    colName = varNames{dataIdx(c)};

    aMask = strcmp(string(anchorTable.Column), string(colName));
    if any(aMask)
        aRow = anchorTable(find(aMask, 1, 'first'), :);
        anchorOK = logical(aRow.AnchorOK);
        P0_h = aRow.Period_hours;
        LightDur_h = aRow.LightDur_h;
        anchorModeUsed = string(aRow.AnchorModeUsed);
    else
        anchorOK = false;
        P0_h = NaN;
        LightDur_h = NaN;
        anchorModeUsed = "NA";
    end

    chosenMinKey = 'Min360';
    chosenMinMinutes = 360;
    crossReason = '';
    crossMetrics = '';

    if anchorOK
        cMask = strcmp(string(T_cross.Column), string(colName));
        if any(cMask)
            cRow = T_cross(find(cMask, 1, 'first'), :);
            nRMS_cross = cRow.nRMS_60vs360;
            dVar_cross = cRow.DeltaVarExpl_60minus360;
            KL360 = cRow.KLikelyMax_Min360;
            KL60  = cRow.KLikelyMax_Min60;

            if (nRMS_cross >= nRMS_high) || (dVar_cross >= dVar_high)
                chosenMinKey = 'Min60';
                chosenMinMinutes = 60;
                crossReason = 'Min 60 chosen: cross-min nRMS and/or DeltaVarExpl above high thresholds';
            elseif (nRMS_cross >= nRMS_low) || (dVar_cross >= dVar_low)
                chosenMinKey = 'Min360';
                chosenMinMinutes = 360;
                crossReason = 'Min 360 chosen: cross-min differences moderate (below high thresholds)';
            else
                chosenMinKey = 'Min360';
                chosenMinMinutes = 360;
                crossReason = 'Min 360 chosen: cross-min differences below low thresholds';
            end
            crossMetrics = sprintf('CrossMin: nRMS=%.3f, dVar=%.3f | KLikelyMax(P360)=%.3g, KLikelyMax(P60)=%.3g', ...
                nRMS_cross, dVar_cross, KL360, KL60);
        else
            chosenMinKey = 'Min360';
            chosenMinMinutes = 360;
            crossReason = 'Min 360 chosen: no cross-min QC row';
            crossMetrics = 'CrossMin: NA';
        end
    else
        crossReason = 'Anchor rejected: pass-through';
        crossMetrics = 'CrossMin: NA';
    end

    methodCode = NC.METHOD_SEL;
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
                methodCode = NC.METHOD_SEL;
                withinReason = 'Selective chosen: within-min nRMS and DeltaVarExpl below low thresholds';
            elseif (nRMS >= nRMS_high) || (dVar >= dVar_high)
                methodCode = NC.METHOD_FULL;
                withinReason = 'Full ladder chosen: within-min nRMS and/or DeltaVarExpl above high thresholds';
            else
                methodCode = NC.METHOD_SEL;
                withinReason = 'Selective chosen: within-min disagreement moderate';
            end
            withinMetrics = sprintf('WithinMin(%s): nRMS=%.3f, dVar=%.3f', chosenMinKey, nRMS, dVar);
        else
            methodCode = NC.METHOD_SEL;
            withinReason = sprintf('Selective chosen: no within-min QC row (%s)', chosenMinKey);
            withinMetrics = sprintf('WithinMin(%s): NA', chosenMinKey);
        end
    else
        methodCode = NC.METHOD_SEL;
        withinReason = 'Anchor rejected: default label only';
        withinMetrics = 'WithinMin: NA';
    end

    recResidual = hs_label(methodCode, chosenMinMinutes, NC);

    workbook = fullfile(NC.TS_DIR, NC.TS_RES_DIR, hs_ts_workbook_name(recResidual, fileStem));
    sigHeader = hs_ts_header(recResidual, colName, false, NC);

    reason = strjoin({crossReason, withinReason, crossMetrics, withinMetrics}, ' | ');

    recRows(end+1,:) = {LightDur_h, colName, anchorOK, P0_h, char(anchorModeUsed), recResidual, workbook, sigHeader, reason};
end

recommendationTable = cell2table(recRows, 'VariableNames', recHeader);

%% -------------------- BUILD VALIDATION MANIFEST SHEET -------------------

% v12: this table separates the algorithmic residual recommendation from the
% downstream validation design.  The recommended residual is retained as QC
% metadata, while SEL_P360 is pre-specified as the primary HSub validation
% residual for Raw-vs-HSub ultradian-period validation.
validationManifest = build_hsub_validation_manifest( ...
    recommendationTable, anchorTable, T_cross, T_qc, fileStem, NC, ...
    nRMS_low, nRMS_high, dVar_low, dVar_high);

validationManifestMAT = fullfile(reportsFolder, 'HS_ValidationManifest.mat');
try
    save(validationManifestMAT, 'validationManifest');
catch ME
    fprintf('Validation manifest MAT save failed: %s\n', ME.message);
end

%% ---------------------------- README + INDEX ----------------------------

readmeLines = {
    sprintf('%s summary', PIPELINE_NAME)
    sprintf('Version: %s', PIPELINE_VERSION)
    sprintf('Input file: %s', inputBaseWithExt)
    ''
    'Overview'
    '- Anchor detection scans 22 to 28 h with block-shuffle p-value gating'
    '- Full ladder removes k=1..K down to Min 360 min and Min 60 min'
    '- Selective removes k that behave like harmonics across windows (plus k=1)'
    '- Anchor basis auto-routes fixed vs time-varying using windowed period diagnostics'
    '- Selective includes locked vs free evidence guard (median dR2_freeMinusLocked <= 0.02)'
    ''
    'Outputs'
    sprintf('- Time series: %s/%s (residual) and %s/%s (removed)', NC.TS_DIR, NC.TS_RES_DIR, NC.TS_DIR, NC.TS_REM_DIR)
    sprintf('- Scalograms (if enabled): %s/%s and %s/%s', NC.SCALO_DIR, NC.SCALO_RES_DIR, NC.SCALO_DIR, NC.SCALO_REM_DIR)
    sprintf('- Anchor figures (if enabled): %s', NC.FIG_ANCHOR_DIR)
    ''
    'Residual labels'
    '- FL_P360, FL_P60, SEL_P360, SEL_P60'
    ''
    'v12 validation design'
    '- The Recommend sheet still records the algorithmic residual recommendation'
    '- The ValidationManifest sheet pre-specifies SEL_P360 as the primary HSub validation residual'
    '- SEL_P60 is retained as secondary validation; FL_P360 and FL_P60 are retained as aggressive sensitivity outputs'
    '- Full-ladder recommendations are retained as QC/sensitivity flags, not as the default biological validation mode'
    ''
    'QC interpretation'
    '- AnchorOK means circadian-scale support (pBlock<0.05 plus viability guards)'
    '- AnchorModeUsed is fixed vs time-varying anchor phase basis'
    '- P0w_* summarise windowed anchor period variability used for auto routing'
    ''
    'Selective vs full ladder QC (within-min)'
    sprintf('- nRMS: <%.2f low, %.2f to %.2f moderate, >=%.2f high', nRMS_low, nRMS_low, nRMS_high, nRMS_high)
    sprintf('- DeltaVarExpl: <%.2f small, %.2f to %.2f moderate, >=%.2f large', dVar_low, dVar_low, dVar_high, dVar_high)
    ''
    'Cross-min QC (full ladder: Min 60 vs Min 360)'
    '- Large nRMS_60vs360 and DeltaVarExpl_60minus360 implies higher-order harmonic contribution'
    ''
    'Recommend sheet'
    '- Picks one residual per column using cross-min QC and within-min QC'
    '- Workbook and SignalColumnHeader match TS/RES naming'
    '- This is an algorithmic QC recommendation, not the default downstream ultradian validation rule'
    ''
    'ValidationManifest sheet'
    '- PrimaryValidationResidual is SEL_P360 by design for the Raw-vs-HSub validation pipeline'
    '- RecommendedResidual is retained for QC, sensitivity interpretation, and auditability'
    '- FullLadderRecommendedFlag indicates cases where aggressive subtraction was algorithmically preferred'
    ''
    'Excluded sheet'
    '- Data columns excluded in the optional preflight exclude step'
    ''
    'Susceptibility (detail workbook)'
    '- AttenuationRatio = GWS_Min60 / GWS_Min360'
    '- <0.30 harmonic-susceptible, 0.30 to 0.70 ambiguous, >0.70 harmonic-robust'
    ''
    'Errors'
    '- Caught per-column exceptions during processing'
    };

indexRows = {
    'Sheet','Contents'
    SH.SUM_README,'How to interpret outputs and QC'
    SH.SUM_INDEX,'Sheet map'
    SH.SUM_EXCL,'User-excluded columns'
    SH.SUM_GAP,'Missingness and interpolation actions'
    SH.SUM_ANCH,'Anchor scan, gating, and auto anchor mode'
    SH.SUM_FL360,'Full ladder summary (Min 360 min)'
    SH.SUM_FL60,'Full ladder summary (Min 60 min)'
    SH.SUM_SEL360,'Selective summary (Min 360 min)'
    SH.SUM_SEL60,'Selective summary (Min 60 min)'
    SH.SUM_QC_SF_360,'Selective vs full ladder QC (Min 360 min)'
    SH.SUM_QC_SF_60,'Selective vs full ladder QC (Min 60 min)'
    SH.SUM_QC_CM,'Full ladder cross-min QC (Min 60 vs Min 360)'
    SH.SUM_REC,'Recommended residual per column'
    SH.SUM_VALMAN,'v12 downstream validation manifest: SEL_P360 primary validation plus QC/sensitivity metadata'
    SH.SUM_ERR,'Caught processing errors'
    };

%% ------------------------ WRITE SUMMARY WORKBOOK ------------------------

enforce_output_path_length(summaryXLSX, 'Summary workbook');

writecell(readmeLines(:), summaryXLSX, 'Sheet', SH.SUM_README);
writecell(indexRows,      summaryXLSX, 'Sheet', SH.SUM_INDEX);

writetable(excludedTable, summaryXLSX, 'Sheet', SH.SUM_EXCL);
writetable(gapTable,      summaryXLSX, 'Sheet', SH.SUM_GAP);
writetable(anchorTable,   summaryXLSX, 'Sheet', SH.SUM_ANCH);

writetable(cell2table(fullReport.Min360, 'VariableNames', fullHeader), summaryXLSX, 'Sheet', SH.SUM_FL360);
writetable(cell2table(fullReport.Min60,  'VariableNames', fullHeader), summaryXLSX, 'Sheet', SH.SUM_FL60);
writetable(cell2table(selReport.Min360,  'VariableNames', selHeader),  summaryXLSX, 'Sheet', SH.SUM_SEL360);
writetable(cell2table(selReport.Min60,   'VariableNames', selHeader),  summaryXLSX, 'Sheet', SH.SUM_SEL60);
writetable(cell2table(qcSelFull.Min360,  'VariableNames', qcSelFullHeader), summaryXLSX, 'Sheet', SH.SUM_QC_SF_360);
writetable(cell2table(qcSelFull.Min60,   'VariableNames', qcSelFullHeader), summaryXLSX, 'Sheet', SH.SUM_QC_SF_60);
writetable(T_cross, summaryXLSX, 'Sheet', SH.SUM_QC_CM);

writetable(recommendationTable, summaryXLSX, 'Sheet', SH.SUM_REC);
writetable(validationManifest, summaryXLSX, 'Sheet', SH.SUM_VALMAN);
writetable(errorsTable, summaryXLSX, 'Sheet', SH.SUM_ERR);

%% ------------------------- WRITE DETAIL WORKBOOK ------------------------

enforce_output_path_length(detailXLSX, 'Detail workbook');

detailIndex = {
    'Sheet','Contents'
    SH.DET_INDEX,'Sheet map'
    SH.DET_SELW360,'Selective per-window fits (Min 360 min; locked vs free evidence)'
    SH.DET_SELW60,'Selective per-window fits (Min 60 min; locked vs free evidence)'
    SH.DET_SUSC,'Global wavelet susceptibility (Min 360 vs Min 60)'
    };
writecell(detailIndex, detailXLSX, 'Sheet', SH.DET_INDEX);
writetable(cell2table(selWinReport.Min360, 'VariableNames', selWinHeader), detailXLSX, 'Sheet', SH.DET_SELW360);
writetable(cell2table(selWinReport.Min60,  'VariableNames', selWinHeader), detailXLSX, 'Sheet', SH.DET_SELW60);
writetable(T_sus, detailXLSX, 'Sheet', SH.DET_SUSC);

fprintf('\nDone.\n');
fprintf('Reports:\n  %s\n  %s\n', summaryXLSX, detailXLSX);
end


%% ------------------------------------------------------------------------
% v12 validation-manifest helpers
% ------------------------------------------------------------------------
function validationManifest = build_hsub_validation_manifest(recommendationTable, anchorTable, T_cross, T_qc, fileStem, NC, nRMS_low, nRMS_high, dVar_low, dVar_high)

    manifestHeader = { ...
        'FileStem', ...
        'LightDur_h', ...
        'Column', ...
        'AnchorOK', ...
        'AnchorPeriod_h', ...
        'AnchorModeUsed', ...
        'RecommendedResidual', ...
        'RecommendedWorkbook', ...
        'RecommendedSignalColumnHeader', ...
        'RecommendationReason', ...
        'PrimaryValidationResidual', ...
        'PrimaryValidationWorkbook', ...
        'PrimaryValidationSignalColumnHeader', ...
        'SecondaryValidationResidual', ...
        'SecondaryValidationWorkbook', ...
        'SecondaryValidationSignalColumnHeader', ...
        'SensitivityResiduals', ...
        'FullLadderRecommendedFlag', ...
        'FullLadderRecommendedReason', ...
        'CrossMinRecommendation', ...
        'CrossMinRisk', ...
        'WithinP360Risk', ...
        'WithinP60Risk', ...
        'ValidationUse', ...
        'HSubQCFlag', ...
        'QCNotes'};

    if isempty(recommendationTable) || height(recommendationTable) == 0
        validationManifest = cell2table(cell(0, numel(manifestHeader)), 'VariableNames', manifestHeader);
        return;
    end

    primaryResidual   = hs_label(NC.METHOD_SEL, 360, NC);  % SEL_P360
    secondaryResidual = hs_label(NC.METHOD_SEL, 60,  NC);  % SEL_P60
    sensitivityResiduals = strjoin({hs_label(NC.METHOD_FULL, 360, NC), hs_label(NC.METHOD_FULL, 60, NC)}, '; ');

    rows = cell(height(recommendationTable), numel(manifestHeader));

    for ii = 1:height(recommendationTable)
        colName = string(table_get_scalar(recommendationTable, ii, 'Column', ""));

        LightDur_h = table_get_scalar(recommendationTable, ii, 'LightDur_h', NaN);
        anchorOK   = logical(table_get_scalar(recommendationTable, ii, 'AnchorOK', false));
        P0_h       = table_get_scalar(recommendationTable, ii, 'Period_hours', NaN);
        anchorMode = string(table_get_scalar(recommendationTable, ii, 'AnchorModeUsed', "NA"));

        recResidual = string(table_get_scalar(recommendationTable, ii, 'RecommendedResidual', "NA"));
        recWorkbook = string(table_get_scalar(recommendationTable, ii, 'Workbook', ""));
        recHeader   = string(table_get_scalar(recommendationTable, ii, 'SignalColumnHeader', ""));
        recReason   = string(table_get_scalar(recommendationTable, ii, 'Reason', ""));

        % Prefer anchor-table values if present; these are the canonical anchor statistics.
        if ~isempty(anchorTable) && height(anchorTable) > 0
            aIdx = find(strcmp(string(anchorTable.Column), colName), 1, 'first');
            if ~isempty(aIdx)
                LightDur_h = table_get_scalar(anchorTable, aIdx, 'LightDur_h', LightDur_h);
                anchorOK   = logical(table_get_scalar(anchorTable, aIdx, 'AnchorOK', anchorOK));
                P0_h       = table_get_scalar(anchorTable, aIdx, 'Period_hours', P0_h);
                anchorMode = string(table_get_scalar(anchorTable, aIdx, 'AnchorModeUsed', anchorMode));
            end
        end

        primaryWorkbook = string(fullfile(NC.TS_DIR, NC.TS_RES_DIR, hs_ts_workbook_name(primaryResidual, fileStem)));
        primaryHeader   = string(hs_ts_header(primaryResidual, char(colName), false, NC));
        secondaryWorkbook = string(fullfile(NC.TS_DIR, NC.TS_RES_DIR, hs_ts_workbook_name(secondaryResidual, fileStem)));
        secondaryHeader   = string(hs_ts_header(secondaryResidual, char(colName), false, NC));

        [crossRec, crossRisk, crossNote] = validation_cross_min_fields(T_cross, colName, nRMS_low, nRMS_high, dVar_low, dVar_high);
        [risk360, qNote360] = validation_within_min_risk(T_qc, 'Min360', colName);
        [risk60,  qNote60]  = validation_within_min_risk(T_qc, 'Min60',  colName);

        fullFlag = startsWith(recResidual, string(NC.METHOD_FULL));
        fullReason = "";
        if fullFlag
            fullReason = "Algorithmic Recommend sheet selected Full Ladder because within-min Selective-vs-Full disagreement crossed the high threshold.";
        else
            fullReason = "Algorithmic Recommend sheet did not select Full Ladder.";
        end

        if ~anchorOK
            validationUse = "No supportive HSub validation: anchor rejected or unavailable; residual is pass-through for this column.";
            qcFlag = "AnchorRejected";
        else
            validationUse = "Use SEL_P360 as primary HSub validation residual; retain Recommend/Full Ladder fields as QC and sensitivity metadata.";
            qcFlag = classify_manifest_qc(fullFlag, crossRisk, risk360, risk60);
        end

        qcNotes = strjoin(string({ ...
            char(crossNote), ...
            char(qNote360), ...
            char(qNote60), ...
            sprintf('Thresholds: within/cross nRMS low=%.3g high=%.3g; dVar low=%.3g high=%.3g', nRMS_low, nRMS_high, dVar_low, dVar_high)}), ' | ');

        rows(ii,:) = { ...
            fileStem, ...
            LightDur_h, ...
            char(colName), ...
            anchorOK, ...
            P0_h, ...
            char(anchorMode), ...
            char(recResidual), ...
            char(recWorkbook), ...
            char(recHeader), ...
            char(recReason), ...
            char(primaryResidual), ...
            char(primaryWorkbook), ...
            char(primaryHeader), ...
            char(secondaryResidual), ...
            char(secondaryWorkbook), ...
            char(secondaryHeader), ...
            sensitivityResiduals, ...
            fullFlag, ...
            char(fullReason), ...
            char(crossRec), ...
            char(crossRisk), ...
            char(risk360), ...
            char(risk60), ...
            char(validationUse), ...
            char(qcFlag), ...
            char(qcNotes)};
    end

    validationManifest = cell2table(rows, 'VariableNames', manifestHeader);
end

function [crossRec, crossRisk, crossNote] = validation_cross_min_fields(T_cross, colName, nRMS_low, nRMS_high, dVar_low, dVar_high)
    crossRec = "NA";
    crossRisk = "NA";
    crossNote = "CrossMin: no row";

    if isempty(T_cross) || height(T_cross) == 0
        return;
    end

    idx = find(strcmp(string(T_cross.Column), string(colName)), 1, 'first');
    if isempty(idx)
        return;
    end

    crossRec = string(table_get_scalar(T_cross, idx, 'Recommendation', "NA"));
    nRMS = table_get_scalar(T_cross, idx, 'nRMS_60vs360', NaN);
    dVar = table_get_scalar(T_cross, idx, 'DeltaVarExpl_60minus360', NaN);

    if isfinite(nRMS) || isfinite(dVar)
        if (isfinite(nRMS) && nRMS >= nRMS_high) || (isfinite(dVar) && dVar >= dVar_high)
            crossRisk = "High";
        elseif (isfinite(nRMS) && nRMS >= nRMS_low) || (isfinite(dVar) && dVar >= dVar_low)
            crossRisk = "Moderate";
        else
            crossRisk = "Low";
        end
        crossNote = sprintf('CrossMin: Recommendation=%s, nRMS=%.4g, dVar=%.4g', char(crossRec), nRMS, dVar);
    else
        crossRisk = "NA";
        crossNote = sprintf('CrossMin: Recommendation=%s, metrics unavailable', char(crossRec));
    end
end

function [risk, note] = validation_within_min_risk(T_qc, minKey, colName)
    risk = "NA";
    note = sprintf('Within%s: no row', minKey);

    if isempty(T_qc) || ~isfield(T_qc, minKey)
        return;
    end

    qTbl = T_qc.(minKey);
    if isempty(qTbl) || height(qTbl) == 0
        return;
    end

    idx = find(strcmp(string(qTbl.Column), string(colName)), 1, 'first');
    if isempty(idx)
        return;
    end

    risk = string(table_get_scalar(qTbl, idx, 'Risk', "NA"));
    nRMS = table_get_scalar(qTbl, idx, 'nRMS', NaN);
    dVar = table_get_scalar(qTbl, idx, 'DeltaVarExpl_FullMinusSel', NaN);
    note = sprintf('Within%s: Risk=%s, nRMS=%.4g, dVar=%.4g', minKey, char(risk), nRMS, dVar);
end

function qcFlag = classify_manifest_qc(fullFlag, crossRisk, risk360, risk60)
    risks = upper(string({crossRisk, risk360, risk60}));
    if fullFlag || any(risks == "HIGH")
        qcFlag = "High";
    elseif any(risks == "MODERATE") || any(risks == "REVIEW")
        qcFlag = "Moderate";
    elseif any(risks == "LOW")
        qcFlag = "Low";
    else
        qcFlag = "Unclassified";
    end
end

function value = table_get_scalar(T, rowIdx, varName, defaultValue)
    value = defaultValue;
    if isempty(T) || height(T) < rowIdx || ~ismember(varName, T.Properties.VariableNames)
        return;
    end

    try
        raw = T{rowIdx, varName};
    catch
        return;
    end

    if iscell(raw)
        if isempty(raw) || isempty(raw{1})
            return;
        end
        value = raw{1};
    elseif ischar(raw)
        if isempty(raw)
            return;
        end
        value = raw;
    elseif isstring(raw)
        if isempty(raw)
            return;
        end
        value = raw(1);
    elseif iscategorical(raw)
        if isempty(raw)
            return;
        end
        value = string(raw(1));
    elseif numel(raw) >= 1
        value = raw(1);
    else
        value = raw;
    end
end

%% ------------------------------------------------------------------------
% Preflight + UI helpers
% ------------------------------------------------------------------------
function [fileList, inputParent] = select_input_files(runMode)
    fileList = {};
    inputParent = '';

    if strcmpi(runMode, 'Single file')
        [fileName, filePath] = uigetfile({'*.xlsx','Excel files (*.xlsx)'}, 'Select input .xlsx file');
        if isequal(fileName, 0)
            return;
        end
        fileList = {fullfile(filePath, fileName)};
        inputParent = filePath;
    else
        [files, filePath] = uigetfile({'*.xlsx','Excel files (*.xlsx)'}, ...
            'Select input .xlsx file(s)', ...
            'MultiSelect', 'on');
        if isequal(files, 0)
            return;
        end
        if ischar(files)
            files = {files};
        end
        files = files(:);
        fileList = cell(numel(files),1);
        for i = 1:numel(files)
            fileList{i} = fullfile(filePath, files{i});
        end
        inputParent = filePath;
    end
end

function thresholds = ask_selective_thresholds_dialog()
    prompt = { ...
        'ratioTol (default 0.08): MAD tolerance on harmonic period ratios; smaller is stricter', ...
        'Rcrit (default 0.55): phase-locking threshold; larger requires stronger locking', ...
        'rhoCrit (default 0.45): Spearman coupling threshold (long regime); larger is stricter', ...
        'ampFracCritShort (default 0.04): short regime amplitude fraction; larger is stricter' ...
        };

    dlgTitle = 'Selective thresholds';
    defAns  = {'0.08','0.55','0.45','0.04'};
    answer  = inputdlg(prompt, dlgTitle, [1 95], defAns);

    if isempty(answer)
        thresholds.ratioTol = 0.08;
        thresholds.Rcrit = 0.55;
        thresholds.rhoCrit = 0.45;
        thresholds.ampFracCritShort = 0.04;
        return;
    end

    v = [str2double(answer{1}), str2double(answer{2}), str2double(answer{3}), str2double(answer{4})];
    d = [0.08 0.55 0.45 0.04];
    bad = ~isfinite(v) | v < 0;
    v(bad) = d(bad);

    thresholds.ratioTol = v(1);
    thresholds.Rcrit = v(2);
    thresholds.rhoCrit = v(3);
    thresholds.ampFracCritShort = v(4);
end

function mapPF = preflight_column_mapping_dialog(dataTable, fileIndex, nFiles)

    varNames = dataTable.Properties.VariableNames;
    nVars = numel(varNames);

    [listStr, exactListKeys] = make_preflight_list_labels(varNames);

    timeGuess  = guess_time_column(varNames);
    lightGuess = guess_light_column(varNames);

    if isempty(timeGuess),  timeGuess = 1;      end
    if isempty(lightGuess), lightGuess = nVars; end

    timeGuess  = max(1, min(nVars, timeGuess));
    lightGuess = max(1, min(nVars, lightGuess));

    titleBar = sprintf('Column setup | File %d/%d', fileIndex, nFiles);

    promptTime = {sprintf('File %d/%d', fileIndex, nFiles), 'Select Time column'};

    [timeIdx, okTime] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', promptTime, ...
        'SelectionMode', 'single', ...
        'ListString', listStr, ...
        'InitialValue', timeGuess, ...
        'ListSize', [640 420]);

    if ~okTime
        error('No Time column selected. Preflight cancelled.');
    end

    promptLight = {sprintf('File %d/%d', fileIndex, nFiles), 'Select Light duration (h) column'};

    [lightIdx, okLight] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', promptLight, ...
        'SelectionMode', 'single', ...
        'ListString', listStr, ...
        'InitialValue', lightGuess, ...
        'ListSize', [640 420]);

    if ~okLight
        error('No Light duration column selected. Preflight cancelled.');
    end

    defaultData = setdiff(1:nVars, [timeIdx, lightIdx], 'stable');
    promptData = {sprintf('File %d/%d', fileIndex, nFiles), 'Select Data columns'};

    [dataIdx, okData] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', promptData, ...
        'SelectionMode', 'multiple', ...
        'ListString', listStr, ...
        'InitialValue', defaultData, ...
        'ListSize', [640 420]);

    if ~okData || isempty(dataIdx)
        error('No Data columns selected. Preflight cancelled.');
    end

    dataIdx = setdiff(dataIdx, [timeIdx, lightIdx], 'stable');
    if isempty(dataIdx)
        error('After excluding Time and Light duration, no Data columns remain.');
    end

    % Optional exclude step
    selectedDataNames = varNames(dataIdx);
    selectedDataLabels = listStr(dataIdx);

    excludeList = [{'<<None>>'}; selectedDataLabels(:)];

    promptExclude = { ...
        sprintf('File %d/%d', fileIndex, nFiles), ...
        'Optional: select Data columns to exclude (or keep <<None>>)', ...
        };

    [excludeSel, okExclude] = listdlg( ...
        'Name', titleBar, ...
        'PromptString', promptExclude, ...
        'SelectionMode', 'multiple', ...
        'ListString', excludeList, ...
        'InitialValue', 1, ...
        'ListSize', [700 420]);

    excludedNames = {};
    excludedIdxPre = [];
    if okExclude
        excludeSel = setdiff(excludeSel, 1);
        if ~isempty(excludeSel)
            excludePos = excludeSel - 1;
            excludePos = excludePos(excludePos >= 1 & excludePos <= numel(dataIdx));
            if ~isempty(excludePos)
                excludedIdxPre = dataIdx(excludePos);
                excludedNames = selectedDataNames(excludePos);
                dataIdx(excludePos) = [];
            end
        end
    end

    if isempty(dataIdx)
        error('All Data columns were excluded.');
    end

    mapPF = struct();
    mapPF.nColsPreflight = nVars;
    mapPF.varNamesPreflight = varNames(:)';
    mapPF.listLabelsPreflight = listStr(:)';
    mapPF.listExactKeys = exactListKeys(:)';

    mapPF.timeIdx_pre = timeIdx;
    mapPF.lightIdx_pre = lightIdx;
    mapPF.dataIdx_pre = dataIdx(:)';

    mapPF.timeName_pre = varNames{timeIdx};
    mapPF.lightName_pre = varNames{lightIdx};
    mapPF.dataNames_pre = varNames(dataIdx);

    mapPF.timeCanon = canonicalise_header(varNames{timeIdx});
    mapPF.lightCanon = canonicalise_header(varNames{lightIdx});
    mapPF.dataCanon = cellfun(@canonicalise_header, varNames(dataIdx), 'UniformOutput', false);

    mapPF.excludedIdx_pre = excludedIdxPre(:)';
    mapPF.excludedNames = excludedNames(:)';

    mapPF.fileIndex = fileIndex;
    mapPF.fileCount = nFiles;
end

function outFolder = unique_subfolder_for_file(rootFolder, fileIndex, fileStem)
    safeStem = sanitise_filename(fileStem);
    baseName = sprintf('File_%02d_%s', fileIndex, safeStem);

    rootFolder = char(string(rootFolder));
    while ~isempty(rootFolder) && any(rootFolder(end) == ['/' '\'])
        rootFolder(end) = [];
    end

    [~, rootLeaf] = fileparts(rootFolder);
    if strcmpi(rootLeaf, baseName)
        outFolder = rootFolder;
        return;
    end

    outFolder = fullfile(rootFolder, baseName);

    if ~exist(outFolder, 'dir')
        return;
    end

    k = 2;
    while true
        outFolderCandidate = fullfile(rootFolder, sprintf('%s_%02d', baseName, k));
        if ~exist(outFolderCandidate, 'dir')
            outFolder = outFolderCandidate;
            return;
        end
        k = k + 1;
    end
end

function [timeIdx, lightIdx, dataIdx, excludedNames, notes] = resolve_preflight_mapping(dataTable, mapPF)

    varNamesNow = dataTable.Properties.VariableNames;
    notes = {};

    [timeIdx, noteT] = resolve_one_column(varNamesNow, mapPF.timeName_pre, mapPF.timeCanon, mapPF.timeIdx_pre, 'Time');
    if ~isempty(noteT), notes{end+1} = noteT; end %#ok<AGROW>

    [lightIdx, noteL] = resolve_one_column(varNamesNow, mapPF.lightName_pre, mapPF.lightCanon, mapPF.lightIdx_pre, 'Light duration (h)');
    if ~isempty(noteL), notes{end+1} = noteL; end %#ok<AGROW>

    dataIdx = [];
    missingData = {};
    dataNames_pre = mapPF.dataNames_pre;
    dataCanon_pre = mapPF.dataCanon;
    dataIdx_pre = mapPF.dataIdx_pre;

    nData = numel(dataNames_pre);
    for i = 1:nData
        thisName = dataNames_pre{i};
        thisCanon = dataCanon_pre{i};
        thisIdxPre = NaN;
        if i <= numel(dataIdx_pre)
            thisIdxPre = dataIdx_pre(i);
        end

        [idxNow, noteNow] = resolve_one_column(varNamesNow, thisName, thisCanon, thisIdxPre, sprintf('Data[%d]', i));
        if isempty(idxNow)
            missingData{end+1} = thisName; %#ok<AGROW>
        else
            dataIdx(end+1) = idxNow; %#ok<AGROW>
            if ~isempty(noteNow)
                notes{end+1} = noteNow; %#ok<AGROW>
            end
        end
    end

    dataIdx = unique(dataIdx, 'stable');
    dataIdx = setdiff(dataIdx, [timeIdx, lightIdx], 'stable');

    if ~isempty(missingData)
        notes{end+1} = sprintf('Preflight Data columns not found (skipped): %s', strjoin(missingData, ', '));
    end

    if isempty(timeIdx)
        error('Preflight Time column not found: "%s"', mapPF.timeName_pre);
    end
    if isempty(lightIdx)
        error('Preflight Light duration column not found: "%s"', mapPF.lightName_pre);
    end
    if isempty(dataIdx)
        error('No preflight Data columns could be resolved in the current file.');
    end

    if isfield(mapPF, 'excludedNames') && ~isempty(mapPF.excludedNames)
        excludedNames = mapPF.excludedNames;
    else
        excludedNames = {};
    end
end

function [idx, note] = resolve_one_column(varNamesNow, targetName_pre, targetCanon, idxPre, roleLabel)
    idx = [];
    note = '';

    exactMatch = find(strcmp(varNamesNow, targetName_pre), 1, 'first');
    if ~isempty(exactMatch)
        idx = exactMatch;
        return;
    end

    canonNow = cellfun(@canonicalise_header, varNamesNow, 'UniformOutput', false);
    canonMatch = find(strcmp(canonNow, targetCanon), 1, 'first');
    if ~isempty(canonMatch)
        idx = canonMatch;
        note = sprintf('%s resolved by canonical match: "%s" -> "%s"', ...
            roleLabel, targetName_pre, varNamesNow{canonMatch});
        return;
    end

    if ~isempty(idxPre) && isfinite(idxPre) && idxPre >= 1 && idxPre <= numel(varNamesNow)
        idx = idxPre;
        note = sprintf('%s resolved by stored index (%d): "%s"', roleLabel, idxPre, varNamesNow{idxPre});
        return;
    end
end

function [labels, exactKeys] = make_preflight_list_labels(varNames)
    n = numel(varNames);
    labels = cell(n,1);
    exactKeys = cell(n,1);

    rawLabels = cell(n,1);
    for i = 1:n
        rawLabels{i} = char(string(varNames{i}));
        if isempty(strtrim(rawLabels{i}))
            rawLabels{i} = sprintf('(blank_header_col_%d)', i);
        end
    end

    [u, ~, ic] = unique(rawLabels, 'stable');
    counts = accumarray(ic, 1);

    dupCounter = zeros(size(u));
    for i = 1:n
        label = rawLabels{i};
        g = ic(i);
        if counts(g) > 1
            dupCounter(g) = dupCounter(g) + 1;
            label = sprintf('%s  [#%d]', label, dupCounter(g));
        end
        labels{i} = label;
        exactKeys{i} = rawLabels{i};
    end
end

function idx = guess_time_column(varNames)
    idx = [];
    canons = cellfun(@canonicalise_header, varNames, 'UniformOutput', false);

    strong = find(cellfun(@(s) contains(s,'time') && (contains(s,'hr') || contains(s,'hour') || contains(s,'min') || contains(s,'date')), canons), 1, 'first');
    if ~isempty(strong), idx = strong; return; end

    anyTime = find(cellfun(@(s) contains(s,'time') || contains(s,'datetime') || contains(s,'date'), canons), 1, 'first');
    if ~isempty(anyTime), idx = anyTime; end
end

function idx = guess_light_column(varNames)
    idx = [];
    canons = cellfun(@canonicalise_header, varNames, 'UniformOutput', false);

    strong = find(cellfun(@(s) contains(s,'light') && (contains(s,'duration') || contains(s,'h') || contains(s,'hour')), canons), 1, 'first');
    if ~isempty(strong), idx = strong; return; end

    anyLight = find(cellfun(@(s) contains(s,'light'), canons), 1, 'first');
    if ~isempty(anyLight), idx = anyLight; end
end

function c = canonicalise_header(x)
    c = lower(char(string(x)));
    c = regexprep(c, '\s+', '');
    c = regexprep(c, '[^a-z0-9]', '');
end

%% ------------------------------------------------------------------------
% Output path safety helpers
% ------------------------------------------------------------------------
function warn_if_path_long(pathStr, label)
    if nargin < 2 || isempty(label)
        label = 'Path';
    end

    p = char(string(pathStr));
    nChars = numel(p);

    if ispc && nChars > 180
        fprintf('Warning: %s is long (%d chars). Workbook writing can fail on Windows for deep paths.\n', label, nChars);
        fprintf('  %s\n', p);
    end
end

function enforce_output_path_length(pathStr, label)
    if nargin < 2 || isempty(label)
        label = 'Output path';
    end

    p = char(string(pathStr));
    nChars = numel(p);

    if ispc && nChars > 240
        error(['%s path too long for reliable workbook writing on Windows (%d chars): %s\n' ...
               'Choose a shorter output root (for example C:\MATLAB_Out\Run1) and avoid using an existing File_XX_* folder as the root.'], ...
               label, nChars, p);
    end
end

%% ------------------------------------------------------------------------
% Import helpers
% ------------------------------------------------------------------------
function [dataTable, info] = read_input_table_preserve_robust(inputFile)

    if nargin < 1 || isempty(inputFile)
        error('read_input_table_preserve_robust:NoFile', 'Input file path is required.');
    end
    if ~isfile(inputFile)
        error('read_input_table_preserve_robust:FileNotFound', 'File not found: %s', inputFile);
    end

    info = struct();
    info.File = inputFile;
    info.Method = '';
    info.Sheet = '';
    info.OriginalVariableNames = {};

    try
        opts = detectImportOptions(inputFile, 'FileType', 'spreadsheet');
        if isprop(opts, 'VariableNamingRule')
            opts.VariableNamingRule = 'preserve';
        end

        try
            opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA','N/A','',' '});
        catch
        end

        try
            if isprop(opts, 'ImportErrorRule')
                opts.ImportErrorRule = 'omitrow';
            end
        catch
        end
        try
            if isprop(opts, 'MissingRule')
                opts.MissingRule = 'omitrow';
            end
        catch
        end

        dataTable = readtable(inputFile, opts);
        info.Method = 'detectImportOptions+readtable';
        try
            info.Sheet = opts.Sheet;
        catch
            info.Sheet = '';
        end
        info.OriginalVariableNames = dataTable.Properties.VariableNames;
        return;

    catch ME1
        try
            dataTable = readtable(inputFile, 'VariableNamingRule', 'preserve');
            info.Method = 'readtable_preserve_fallback';
            info.Sheet = '';
            info.OriginalVariableNames = dataTable.Properties.VariableNames;
            return;
        catch ME2
            try
                dataTable = readtable(inputFile);
                info.Method = 'readtable_plain_fallback';
                info.Sheet = '';
                info.OriginalVariableNames = dataTable.Properties.VariableNames;
                return;
            catch ME3
                error('read_input_table_preserve_robust:ReadFailed', ...
                    'Failed to read "%s". Primary: %s | Fallback1: %s | Fallback2: %s', ...
                    inputFile, ME1.message, ME2.message, ME3.message);
            end
        end
    end
end

function [T, info] = drop_empty_columns_robust(T)

    info = struct();
    info.DroppedIdx = [];
    info.DroppedNames = {};

    if isempty(T) || width(T) == 0
        return;
    end

    vn = T.Properties.VariableNames;
    dropMask = false(1, width(T));

    for i = 1:width(T)
        x = T{:, i};
        dropMask(i) = is_column_all_missing(x);
    end

    if any(dropMask)
        info.DroppedIdx = find(dropMask);
        info.DroppedNames = vn(dropMask);
        T(:, dropMask) = [];
    end
end

function tfCol = is_column_all_missing(x)
    try
        if isdatetime(x)
            tfCol = all(isnat(x));
            return;
        end
        if isnumeric(x) || islogical(x)
            tfCol = all(isnan(double(x)));
            return;
        end
        if isstring(x)
            tfCol = all(ismissing(x));
            return;
        end
        if iscell(x)
            tfCol = true;
            for k = 1:numel(x)
                v = x{k};
                if isempty(v)
                    continue;
                end
                if ischar(v) || (isstring(v) && isscalar(v))
                    if strlength(string(v)) == 0 || ismissing(string(v))
                        continue;
                    end
                elseif isnumeric(v) && isscalar(v) && isnan(v)
                    continue;
                else
                    tfCol = false;
                    return;
                end
            end
            return;
        end
        sx = string(x);
        tfCol = all(ismissing(sx));
    catch
        tfCol = false;
    end
end

%% ------------------------------------------------------------------------
% Core analysis helpers
% ------------------------------------------------------------------------
function s = tf(x)
    if x, s = 'OK'; else, s = 'NOT OK'; end
end

function v = light_value_repr(lightVec)
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

function [X0, Xh] = build_drift_and_harmonics_design_phi(t_hr, phi0, kList)
    tz = zscore_local(t_hr(:));
    X0 = [ones(size(tz)), tz, tz.^2];

    ph = phi0(:);
    n = numel(ph);
    nK = numel(kList);

    Xh = zeros(n, 2*nK);
    for ii = 1:nK
        kk = kList(ii);
        col = 2*(ii-1) + 1;
        Xh(:, col)   = cos(kk*ph);
        Xh(:, col+1) = sin(kk*ph);
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

function [A, phi, R2] = fit_locked_phase_amp_phase_R2_with_linear(phiVec, t_hr, y)
    ph = phiVec(:);
    t = t_hr(:);
    y = y(:);

    if numel(ph) < 10 || all(~isfinite(y))
        A = NaN; phi = NaN; R2 = NaN;
        return;
    end

    tZ = zscore_local(t);
    X = [cos(ph), sin(ph), ones(size(ph)), tZ];
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
    % prefix is a compact label (e.g. FL_P360, SEL_P60). Folder determines RES vs REM header style.

    NC = hs_naming();

    N = size(dataMatrix, 1);

    % Infer whether this is the removed folder to add REM tag in headers
    [~, leaf] = fileparts(char(string(outFolder)));
    isRemoved = strcmpi(leaf, NC.TS_REM_DIR);

    headers = cell(1, 2 + numel(dataIdx));
    headers{1} = timeName;
    for c = 1:numel(dataIdx)
        colN = varNames{dataIdx(c)};
        headers{1 + c} = hs_ts_header(prefix, colN, isRemoved, NC);
    end
    headers{end} = lightName;

    out = cell(N+1, numel(headers));
    out(1,:) = headers;
    out(2:end,1) = column_to_cell(timeOut);
    out(2:end,2:(1+numel(dataIdx))) = num2cell(dataMatrix);
    out(2:end,end) = column_to_cell(lightOut);

    outFile = fullfile(outFolder, hs_ts_workbook_name(prefix, fileName));
    enforce_output_path_length(outFile, sprintf('TS workbook (%s)', prefix));

    writecell(out, outFile);
    fprintf('Saved TS (%s): %s\n', prefix, outFile);
end

% ---------------------- v10: windowed anchor diagnostics ------------------
function [P0w, tCenter_h, nUsed] = windowed_anchor_periods(t_min, y, winStarts_min, winEnds_min, periods_h, minSamplesPerWindow)
    nW = numel(winStarts_min);
    P0w = NaN(nW,1);
    tCenter_h = NaN(nW,1);
    nUsed = 0;

    for w = 1:nW
        wMask = (t_min >= winStarts_min(w)) & (t_min < winEnds_min(w));
        if sum(wMask) < minSamplesPerWindow
            continue;
        end
        tmw_hr = (t_min(wMask) / 60);
        ymw    = y(wMask);

        [P0_h_w, ~, ~, ~, ~, ~] = regression_anchor_scan(tmw_hr, ymw, periods_h);

        if isfinite(P0_h_w)
            P0w(w) = P0_h_w;
            tCenter_h(w) = 0.5*(winStarts_min(w) + winEnds_min(w))/60;
            nUsed = nUsed + 1;
        end
    end

    ok = isfinite(P0w) & isfinite(tCenter_h);
    P0w = P0w(ok);
    tCenter_h = tCenter_h(ok);
end

function P0w = window_anchor_period_single(t_hr, y, periods_h)
    [P0w, ~, ~, ~, ~, ~] = regression_anchor_scan(t_hr(:), y(:), periods_h);
end

function [phi0, P0_t_h] = build_timevarying_anchor_phase(t_hr, t_min, y, winStarts_min, winEnds_min, periods_h, minSamplesPerWindow, fallbackP0_h)

    [P0w, tCenter_h, nUsed] = windowed_anchor_periods(t_min, y, winStarts_min, winEnds_min, periods_h, minSamplesPerWindow);
    if nUsed < 3
        error('Not enough windowed anchor estimates to build time-varying phase.');
    end

    medP = median(P0w, 'omitnan');
    madP = mad(P0w, 1);
    if ~isfinite(madP) || madP <= 0
        madP = 0.1;
    end
    lo = max(1, medP - 6*madP);
    hi = medP + 6*madP;
    P0wC = min(max(P0w, lo), hi);

    P0_t_h = interp1(tCenter_h, P0wC, t_hr, 'pchip', 'extrap');

    P0_t_h(~isfinite(P0_t_h)) = fallbackP0_h;
    P0_t_h = min(max(P0_t_h, 20), 30);

    omega = 2*pi ./ (P0_t_h + eps);
    phi0 = cumtrapz(t_hr, omega);
    phi0 = phi0 - phi0(1);
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

%% ------------------------------------------------------------------------
% Naming helpers (centralised)
% ------------------------------------------------------------------------
function NC = hs_naming()
    % Method codes
    NC.METHOD_FULL = 'FL';   % Full ladder
    NC.METHOD_SEL  = 'SEL';  % Selective

    % Min-period label codes
    NC.MIN360_LABEL = 'P360';
    NC.MIN60_LABEL  = 'P60';

    % Tags
    NC.RES_TAG = 'RES';
    NC.REM_TAG = 'REM';

    % Folders
    NC.REPORTS_DIR    = 'Reports';
    NC.TS_DIR         = 'TS';
    NC.TS_RES_DIR     = 'RES';
    NC.TS_REM_DIR     = 'REM';

    NC.FIG_ANCHOR_DIR = 'FIG_ANCHOR';

    NC.SCALO_DIR      = 'SCALO';
    NC.SCALO_RES_DIR  = 'RES';
    NC.SCALO_REM_DIR  = 'REM';

    % Workbooks
    NC.SUMMARY_XLSX_NAME = 'HS_Summary.xlsx';
    NC.DETAIL_XLSX_NAME  = 'HS_Detail.xlsx';

    % Figure prefixes
    NC.ANCH_PREFIX  = 'ANCH';
    NC.SCALO_PREFIX = 'SCALO';
end

function SH = hs_sheets()
    % Summary workbook
    SH.SUM_README = 'README';
    SH.SUM_INDEX  = 'Index';
    SH.SUM_EXCL   = 'Excluded';
    SH.SUM_GAP    = 'GapFill';
    SH.SUM_ANCH   = 'Anchor';
    SH.SUM_FL360  = 'FL_P360';
    SH.SUM_FL60   = 'FL_P60';
    SH.SUM_SEL360 = 'SEL_P360';
    SH.SUM_SEL60  = 'SEL_P60';
    SH.SUM_QC_SF_360 = 'QC_SELvsFL_P360';
    SH.SUM_QC_SF_60  = 'QC_SELvsFL_P60';
    SH.SUM_QC_CM     = 'QC_CrossMin';
    SH.SUM_REC       = 'Recommend';
    SH.SUM_VALMAN    = 'ValidationManifest';
    SH.SUM_ERR       = 'Errors';

    % Detail workbook
    SH.DET_INDEX   = 'Index';
    SH.DET_SELW360 = 'SEL_Win_P360';
    SH.DET_SELW60  = 'SEL_Win_P60';
    SH.DET_SUSC    = 'Suscept_P360vsP60';
end

function label = hs_label(methodCode, minPeriodMinutes, NC)
    label = sprintf('%s_%s', methodCode, hs_min_label(minPeriodMinutes, NC));
end

function minLbl = hs_min_label(minPeriodMinutes, NC)
    if isequal(minPeriodMinutes, 360)
        minLbl = NC.MIN360_LABEL;
    elseif isequal(minPeriodMinutes, 60)
        minLbl = NC.MIN60_LABEL;
    else
        minLbl = sprintf('P%d', round(minPeriodMinutes));
    end
end

function fn = hs_ts_workbook_name(label, fileStem)
    fn = sprintf('%s_%s.xlsx', label, fileStem);
end

function hdr = hs_ts_header(label, colName, isRemoved, NC)
    if isRemoved
        hdr = sprintf('%s_%s_%s', label, NC.REM_TAG, colName);
    else
        hdr = sprintf('%s_%s', label, colName);
    end
end

function fn = hs_scalo_fig_name(label, isRemoved, safeCol, NC)
    tag = NC.RES_TAG;
    if isRemoved
        tag = NC.REM_TAG;
    end
    fn = sprintf('%s_%s_%s_%s.jpg', NC.SCALO_PREFIX, tag, label, safeCol);
end