% -------------------------------------------------------------------------
% LB_02_Wavelet_and_AcrossPhotoperiod.m
% -------------------------------------------------------------------------
% Self-contained code covering:
%   Part A — Behavioural wavelet (per photoperiod file)
%   Part B — Across-photoperiod summaries (cohort level)
%
% Developed for MATLAB R2025b. Helper functions sit at the bottom of this
% file.
%
% Pipeline order (important)
%   1. Run LB_01_HarmonicSubtraction.m first.
%   2. Run Part A of this code on the same activity Excel files.
%      Wavelets are computed on RAW activity, and also on the selective
%      Min360 Removed and Residual series from harmonic subtraction
%      (two-panel scalograms; AnchorOK mice only for HSub averages).
%   3. Run Part B once the cohort Summaries/ folder has all
%      AnalysisSummary__*.mat files. Part B builds RAW (exploratory) and
%      residual-validated (confirmatory) tables plus:
%        co-expression / period comparison, period distributions,
%        stitched RAW scalograms, and stitched Removed|Residual scalograms.
%
% Wavelet settings
%   Morlet CWT ('amor'), period limits 60–1590 min, jet scalograms,
%   top-5 peaks from mean log10 power, bands CR20–28 and UR 1–3 / 3–6 /
%   6–12 / 12–18 h.
% -------------------------------------------------------------------------

clearvars; close all; clc;

fprintf('Wavelet + Across photoperiod\n\n');

%% ========================= WHICH PART TO RUN =============================
% Part A = per-file wavelet. Part B = cohort across-photoperiod summaries.
% "Both" runs A then B in one sitting once the dialogs are answered.

partChoice = questdlg( ...
    ['Part A = wavelet on RAW + HSub Removed/Residual (per photoperiod file).\n' ...
     'Part B = across-photoperiod cohort summaries (needs Summaries/).\n' ...
     'Both = Part A then Part B in one sitting.'], ...
    'Wavelet / Across', ...
    'Part A only', 'Part B only', 'Both', 'Both');
if isempty(partChoice)
    fprintf('Cancelled.\n');
    return;
end

cfg = lb_defaults();   % analysis defaults used by both parts

%% ========================= PART A — WAVELET ==============================
% CWT on RAW (primary), plus Removed|Residual scalograms from LB_01 outputs.

if strcmp(partChoice, 'Part A only') || strcmp(partChoice, 'Both')

    fprintf('\n--- Part A: Wavelet on RAW ---\n');

    runMode = questdlg('Choose input mode:', 'Wavelet', ...
        'Single file', 'Multiple files', 'Cancel', 'Single file');
    if isempty(runMode) || strcmpi(runMode, 'Cancel')
        fprintf('Cancelled Part A.\n');
        if strcmp(partChoice, 'Part A only'), return; end
    else
        [fileList, inPath] = lb_select_input_files(runMode);
        if isempty(fileList)
            fprintf('No input file(s) selected for Part A.\n');
            if strcmp(partChoice, 'Part A only'), return; end
        else
            % Need harmonic-subtraction root so Anchor/Recommendation can be attached
            hsubRoot = uigetdir(inPath, ...
                'Select harmonic subtraction output root (folder containing per-file runs)');
            if isequal(hsubRoot, 0)
                fprintf('No harmonic subtraction root selected. Exiting.\n');
                return;
            end

            outRoot = uigetdir(inPath, ...
                'Select wavelet output root (Wavelet will be created here if needed)');
            if isequal(outRoot, 0)
                fprintf('No wavelet output folder selected. Exiting.\n');
                return;
            end
            if ~contains(lower(outRoot), 'wavelet')
                outRoot = fullfile(outRoot, 'Wavelet');
            end
            lb_ensure_dir(outRoot);

            % Summaries sit beside the Wavelet folder (same parent)
            summaryDir = fullfile(fileparts(outRoot), 'Summaries');
            lb_ensure_dir(summaryDir);

            plotMode = questdlg('Figure export mode:', 'Plot mode', ...
                'development', 'publication', 'development');
            if isempty(plotMode), plotMode = 'development'; end

            % Condition groups: define once on the first file, reuse thereafter
            groupsTemplate = [];
            for iF = 1:numel(fileList)
                fprintf('\n=== Wavelet %d/%d: %s ===\n', iF, numel(fileList), fileList{iF});
                try
                    groupsTemplate = lb_wavelet_run_one_file( ...
                        fileList{iF}, hsubRoot, outRoot, summaryDir, ...
                        cfg, plotMode, groupsTemplate);
                catch ME
                    fprintf(2, 'FAILED wavelet: %s\n%s\n', fileList{iF}, ME.message);
                end
            end

            fprintf('\nPart A finished. Summaries folder:\n  %s\n', summaryDir);
        end
    end
end

%% ========================= PART B — ACROSS PHOTOPERIOD ===================
% Needs every AnalysisSummary__*.mat for the cohort already written by Part A.

if strcmp(partChoice, 'Part B only') || strcmp(partChoice, 'Both')

    fprintf('\n--- Part B: Across photoperiod ---\n');

    summaryDir = uigetdir(pwd, 'Select cohort Summaries/ folder (AnalysisSummary__*.mat)');
    if isequal(summaryDir, 0)
        fprintf('No Summaries folder selected. Exiting Part B.\n');
        return;
    end

    [sel, ok] = listdlg('PromptString', 'Select cohort tag:', ...
        'SelectionMode', 'single', 'ListString', cfg.cohortTags);
    if ~ok
        fprintf('No cohort selected. Exiting Part B.\n');
        return;
    end
    cohortTag = cfg.cohortTags{sel};

    outRoot = uigetdir(summaryDir, ...
        'Select cohort results root (parent of AcrossPhotoperiod_*)');
    if isequal(outRoot, 0)
        fprintf('No output root selected. Exiting Part B.\n');
        return;
    end

    plotMode = questdlg('Figure export mode:', 'Plot mode', ...
        'development', 'publication', 'development');
    if isempty(plotMode), plotMode = 'development'; end

    lb_across_run_cohort(summaryDir, outRoot, cohortTag, plotMode, cfg);
    fprintf('\nPart B finished.\n');
end

fprintf('\nDone.\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function cfg = lb_defaults()
% Analysis defaults for wavelet and across-photoperiod steps.
    cfg.cohortTags = {'C57_LP', 'NR2B_LP', 'NR2B_LD_DD'};
    cfg.samplingMinutes = 10;
    cfg.waveletName = 'amor';              % Morlet analytic wavelet
    cfg.periodLimitsMinutes = [60, 1590]; % 1–26.5 h (Core match)
    cfg.scalogramYTicksHours = 0:4:26;
    cfg.topNPeaks = 5;
    cfg.colormap = 'jet';                  % fixed for all scalograms
    cfg.bands.CR_20_28 = [20, 28];
    cfg.bands.UR_1_3   = [1, 3];
    cfg.bands.UR_3_6   = [3, 6];
    cfg.bands.UR_6_12  = [6, 12];
    cfg.bands.UR_12_18 = [12, 18];
    cfg.bands.coexpression = {'CR_20_28', 'UR_1_3', 'UR_3_6', 'UR_6_12', 'UR_12_18'};
    cfg.summaryPrefix = 'AnalysisSummary__';
    cfg.hsubArm = 'Min360';   % selective residual with 360 min period floor
    cfg.hsubValidation.harmonicToleranceFrac = 0.05;  % ±5% of P0/k
    cfg.hsubValidation.maxHarmonicOrder = 6;
    cfg.hsubValidation.urPeriodRangeHours = [1, 18];
    cfg.hsubValidation.topNPeaks = 5;
    cfg.plot.development = struct('dpi', 96, 'format', 'png');
    cfg.plot.publication = struct('dpi', 600, 'format', 'jpeg');
end

%% --------------------- Part A: wavelet one file -------------------------

function groupsTemplate = lb_wavelet_run_one_file(rawFile, hsubRoot, outRoot, ...
        summaryDir, cfg, plotMode, groupsTemplate)
% Wavelet on RAW for one activity file; write scalograms + summary MAT.

    % --- Read Excel and build time / light vectors ---
    [tbl, varNames] = lb_read_excel(rawFile);
    [timeIdx, lightIdx, mouseIdx] = lb_map_columns(varNames);
    if isempty(timeIdx) || isempty(lightIdx) || isempty(mouseIdx)
        error('Could not resolve Time / Light duration / mouse columns.');
    end

    time_hours = tbl{:, timeIdx};
    if ~isnumeric(time_hours), time_hours = str2double(string(time_hours)); end
    time_min = time_hours * 60;
    time_day = time_min / (60 * 24);   % scalogram x-axis in days
    TsMinutes = median(diff(time_min), 'omitnan');
    if ~isfinite(TsMinutes) || TsMinutes <= 0, TsMinutes = cfg.samplingMinutes; end

    lightVec = tbl{:, lightIdx};
    if ~isnumeric(lightVec), lightVec = str2double(string(lightVec)); end
    LightDur_h = median(lightVec, 'omitnan');
    condChangeIdx = find(diff(lightVec) ~= 0);   % vertical lines on scalograms

    [~, fileStem, ~] = fileparts(rawFile);
    runFolder = fullfile(outRoot, fileStem);
    scaloRoot = fullfile(runFolder, 'Scalograms');
    peakFolder = fullfile(runFolder, 'PeriodPeaks');
    lb_ensure_dir(scaloRoot);
    lb_ensure_dir(peakFolder);
    lb_ensure_dir(summaryDir);

    % Condition groups (e.g. F / M). First file: ask n groups → name → pick columns.
    % Later files reuse the same mouse names from the template.
    if isempty(groupsTemplate)
        groups = lb_group_dialog(varNames, mouseIdx);
        groupsTemplate = groups;
    else
        groups = lb_resolve_groups(varNames, groupsTemplate);
        fprintf('Reusing %d condition group(s).\n', numel(groups));
    end

    % Attach harmonic-subtraction Anchor/Recommendation if present (validation metadata only)
    hsubSummary = fullfile(hsubRoot, fileStem, 'Reports', 'HarmonicRemoval_Summary.xlsx');
    hsubVal = struct('available', isfile(hsubSummary), 'path', hsubSummary);
    if hsubVal.available
        try
            hsubVal.anchor = readtable(hsubSummary, 'Sheet', 'Anchor_Report');
            hsubVal.recommendation = readtable(hsubSummary, 'Sheet', 'Recommendation');
        catch
            hsubVal.available = false;
            warning('Could not read harmonic subtraction summary for %s', fileStem);
        end
    else
        warning('Harmonic subtraction summary not found for %s (hsub.available=false in summary)', fileStem);
    end

    theme = lb_theme(plotMode, cfg);
    % One Morlet filter bank for this recording length / sampling rate
    FB = lb_make_filterbank(height(tbl), TsMinutes, cfg);

    periodAll = {};
    bandAll = {};
    figurePaths = {};

    % --- Group-average scalograms (mean RAW series within each group, then CWT) ---
    for g = 1:numel(groups)
        if isempty(groups(g).colIdx), continue; end
        grpName = groups(g).name;
        grpFolder = fullfile(scaloRoot, ['Group_' lb_sanitise(grpName)]);
        lb_ensure_dir(grpFolder);

        avgSignal = mean(lb_stack_signals(tbl, groups(g).colIdx), 2);
        [wt, periods_hours] = lb_cwt(avgSignal, FB);
        powerSpec = mean(abs(wt).^2, 2);   % time-averaged power spectrum

        outFile = fullfile(grpFolder, sprintf('Scalogram_Average_%s_L%g.%s', ...
            lb_sanitise(grpName), LightDur_h, theme.format));
        lb_plot_scalogram(wt, periods_hours, time_day, condChangeIdx, ...
            sprintf('RAW | Average | %s | n=%d | L=%g h', grpName, numel(groups(g).colIdx), LightDur_h), ...
            outFile, theme, cfg);
        figurePaths{end+1} = outFile; 

        % Top peaks + band powers from the average spectrum
        peaks = lb_find_peaks(powerSpec, periods_hours, cfg.topNPeaks);
        if height(peaks) > 0
            peaks.SignalID = repmat(string(sprintf('Average_%s', lb_sanitise(grpName))), height(peaks), 1);
            peaks.Group = repmat(string(grpName), height(peaks), 1);
            peaks.LightDuration_h = repmat(LightDur_h, height(peaks), 1);
            peaks.N_Mice = repmat(numel(groups(g).colIdx), height(peaks), 1);
            periodAll{end+1} = peaks; 
        end
        bp = lb_band_power(powerSpec, periods_hours, cfg);
        bpRow = struct2table(bp);
        bpRow.SignalID = string(sprintf('Average_%s', lb_sanitise(grpName)));
        bpRow.Group = string(grpName);
        bpRow.LightDuration_h = LightDur_h;
        bpRow.N_Mice = numel(groups(g).colIdx);
        bandAll{end+1} = bpRow; 
    end

    % --- Individual scalograms (same CWT + peak / band summaries per mouse) ---
    for g = 1:numel(groups)
        grpName = groups(g).name;
        indivFolder = fullfile(scaloRoot, ['Group_' lb_sanitise(grpName)], 'Individuals');
        lb_ensure_dir(indivFolder);
        for ii = 1:numel(groups(g).colIdx)
            colIdx = groups(g).colIdx(ii);
            signalID = varNames{colIdx};
            signal = lb_prepare_signal(tbl{:, colIdx});
            [wt, periods_hours] = lb_cwt(signal, FB);
            powerSpec = mean(abs(wt).^2, 2);

            outScalo = fullfile(indivFolder, sprintf('Scalogram_%s_%s.%s', ...
                lb_sanitise(grpName), lb_sanitise(signalID), theme.format));
            lb_plot_scalogram(wt, periods_hours, time_day, condChangeIdx, ...
                sprintf('RAW | %s | %s', signalID, grpName), outScalo, theme, cfg);
            figurePaths{end+1} = outScalo; 

            peaks = lb_find_peaks(powerSpec, periods_hours, cfg.topNPeaks);
            if height(peaks) > 0
                peaks.SignalID = repmat(string(sprintf('%s_%s', grpName, signalID)), height(peaks), 1);
                peaks.Group = repmat(string(grpName), height(peaks), 1);
                peaks.LightDuration_h = repmat(LightDur_h, height(peaks), 1);
                peaks.N_Mice = ones(height(peaks), 1);
                periodAll{end+1} = peaks; 
            end

            bp = lb_band_power(powerSpec, periods_hours, cfg);
            bpRow = struct2table(bp);
            bpRow.SignalID = string(sprintf('%s_%s', grpName, signalID));
            bpRow.Group = string(grpName);
            bpRow.LightDuration_h = LightDur_h;
            bpRow.N_Mice = 1;
            bandAll{end+1} = bpRow; 
        end
    end

    % --- HSub Removed | Residual scalograms (selective Min360; AnchorOK only) ---
    hsubFigPaths = lb_hsub_wavelet_scalograms( ...
        hsubRoot, fileStem, groups, varNames, hsubVal, ...
        time_day, condChangeIdx, LightDur_h, FB, theme, cfg, scaloRoot);
    figurePaths = [figurePaths, hsubFigPaths];

    periodTable = lb_vertcat_tables(periodAll);
    bandTable = lb_vertcat_tables(bandAll);
    writetable(periodTable, fullfile(peakFolder, sprintf('PeriodPeaks_%s.xlsx', fileStem)));

    meta = struct();
    meta.varNames = varNames;
    meta.timeIdx = timeIdx;
    meta.lightIdx = lightIdx;
    meta.mouseIdx = mouseIdx;
    meta.fileStem = fileStem;
    meta.lightDurationHours = LightDur_h;
    meta.cohort = '';

    % Attach mouse column names to groups for residual matching later
    for g = 1:numel(groups)
        groups(g).colNames = varNames(groups(g).colIdx);
    end

    % Summary MAT carries tables + paths into the across-photoperiod step
    payload = struct();
    payload.meta = meta;
    payload.groups = groups;
    payload.hsub = hsubVal;
    payload.periodTable = periodTable;
    payload.bandPower = bandTable;
    payload.figurePaths = figurePaths;
    payload.paths = struct('raw', rawFile, 'hsub', hsubRoot, 'wavelet', runFolder);
    payload.analysisNote = ['Primary wavelet on RAW; also HSub Removed|Residual ' ...
        'scalograms (SEL_P360). Confirmatory UR stats use residual in Part B.'];
    payload.summary_version = '1.0';
    payload.fileStem = fileStem;
    payload.created = datetime('now');

    summaryPath = fullfile(summaryDir, [cfg.summaryPrefix fileStem '.mat']);
    save(summaryPath, '-struct', 'payload', '-v7.3');
    lb_update_summaries_index(summaryDir, fileStem, LightDur_h, summaryPath);
    fprintf('Summary written: %s\n', summaryPath);
    fprintf('Wavelet complete: %s\n', runFolder);
end

%% --------------------- Part B: across photoperiod -----------------------

function lb_across_run_cohort(summaryDir, outputFolder, cohortTag, plotMode, cfg)
% Cohort-level RAW (exploratory) + residual-validated UR (confirmatory).
%
% Tables  — period peaks, co-expression ratios, validation log
% Figures — co-expression lines, period scatter, per-group distributions,
%           stitched RAW scalograms, stitched Removed|Residual scalograms
%
% Photoperiods are ordered by light duration from the summary MAT files.

    outFolder = fullfile(outputFolder, ['AcrossPhotoperiod_' cohortTag]);
    figFolder = fullfile(outFolder, 'Figures');
    tblFolder = fullfile(outFolder, 'Tables');
    lb_ensure_dir(figFolder);
    lb_ensure_dir(tblFolder);

    entries = lb_load_summaries(summaryDir);   % one struct per photoperiod, sorted by L
    theme = lb_theme(plotMode, cfg);

    % --- RAW exploratory layer (from Part A summary tables) ---
    allPeriod = vertcat(entries.periodTable);
    allBand = vertcat(entries.bandPower);
    ratioRAW = lb_coexpression(allBand);

    % --- Validated confirmatory layer (residual CWT; AnchorOK only) ---
    [hsubBand, hsubPeriod, hsubValLog] = lb_hsub_validated_ur(entries, cfg);
    ratioHsub = lb_coexpression(hsubBand);

    writetable(allPeriod, fullfile(tblFolder, 'PeriodPeaks_all_photoperiods_RAW.xlsx'));
    writetable(ratioRAW, fullfile(tblFolder, 'Coexpression_ratios_RAW.xlsx'));
    writetable(hsubBand, fullfile(tblFolder, 'BandPower_hsub_validated.xlsx'));
    writetable(hsubPeriod, fullfile(tblFolder, 'PeriodPeaks_hsub_validated.xlsx'));
    writetable(hsubValLog, fullfile(tblFolder, 'PeriodPeaks_hsub_validation_log.xlsx'));
    writetable(ratioHsub, fullfile(tblFolder, 'Coexpression_ratios_hsub_validated.xlsx'));

    % Figures: RAW first, then confirmatory residual-validated
    lb_plot_coexpression(ratioRAW, fullfile(figFolder, ['Coexpression_CR_UR_RAW.' theme.format]), ...
        'Co-expression RAW (exploratory): CR vs UR band power', 'CR_{20-28} / UR_{total}');
    lb_plot_period_scatter(allPeriod, fullfile(figFolder, ['Period_comparison_RAW.' theme.format]), ...
        'Period comparison RAW (exploratory): rank-1 peak', 'Dominant period (h) — rank 1', [0 26]);

    lb_plot_coexpression(ratioHsub, fullfile(figFolder, ['Coexpression_CR_UR_hsub_validated.' theme.format]), ...
        'Co-expression residual-validated: RAW CR vs residual UR', 'CR_{20-28} (RAW) / UR_{total} (residual)');
    lb_plot_period_scatter(hsubPeriod, fullfile(figFolder, ['Period_comparison_hsub_validated.' theme.format]), ...
        'Residual-validated ultradian periods (confirmatory): rank-1 UR peak', ...
        'Validated UR period (h) — rank 1', [0 18]);

    % Period distributions (one box+jitter figure per condition group)
    % — clearer than a combined scatter when many mice share similar periods
    lb_plot_period_distributions(allPeriod, figFolder, 'Period_distribution_RAW', theme, ...
        'RAW rank-1 period distribution', 'Dominant period (h) — rank 1', [0 26]);
    lb_plot_period_distributions(hsubPeriod, figFolder, 'Period_distribution_hsub_validated', theme, ...
        'Residual-validated UR period distribution', 'Validated UR period (h) — rank 1', [0 18]);

    % Stitched scalograms: concatenate group-mean series across photoperiods
    % in light-duration order, then run one CWT on the long concatenated signal
    lb_plot_stitched_scalograms(entries, figFolder, cohortTag, theme, cfg);
    lb_plot_stitched_hsub_scalograms(entries, figFolder, cohortTag, theme, cfg);

    fprintf('Across-photoperiod complete: %s (%d photoperiods)\n', outFolder, numel(entries));
    fprintf('  Residual-validated UR rows: band=%d, period peaks=%d\n', height(hsubBand), height(hsubPeriod));
end

function entries = lb_load_summaries(summaryDir)
% Load every AnalysisSummary__*.mat and sort by light duration (photoperiod axis).
    mats = dir(fullfile(summaryDir, 'AnalysisSummary__*.mat'));
    if isempty(mats)
        error('No AnalysisSummary__*.mat found in %s', summaryDir);
    end
    entries = struct('fileStem', {}, 'lightHours', {}, 'periodTable', {}, ...
        'bandPower', {}, 'meta', {}, 'groups', {}, 'paths', {}, 'hsub', {});
    for i = 1:numel(mats)
        S = load(fullfile(mats(i).folder, mats(i).name));
        entries(i).fileStem = S.fileStem;
        entries(i).lightHours = S.meta.lightDurationHours;
        entries(i).periodTable = S.periodTable;
        entries(i).bandPower = S.bandPower;
        entries(i).meta = S.meta;
        entries(i).hsub = S.hsub;
        entries(i).groups = S.groups;
        entries(i).paths = S.paths;
    end
    [~, ord] = sort([entries.lightHours]);
    entries = entries(ord);
end

function [bandTable, periodTable, validationLog] = lb_hsub_validated_ur(entries, cfg)
% Confirmatory layer: UR from selective Min360 residual CWT; CR from RAW CWT.
% Only AnchorOK mice. Peaks within ±5% of P0/k (k=1..6) are excluded.

    bandRows = {};
    periodRows = {};
    logRows = {};
    armKey = cfg.hsubArm;   % 'Min360' → Residual_Selective_Min360_*.xlsx

    for i = 1:numel(entries)
        if ~isfield(entries(i), 'hsub') || ~entries(i).hsub.available
            warning('Harmonic subtraction summary missing for %s; skipped.', entries(i).fileStem);
            continue;
        end
        if ~isfield(entries(i).paths, 'hsub')
            warning('paths.hsub missing for %s; skipped.', entries(i).fileStem);
            continue;
        end

        % Selective residual for this photoperiod (360 min period floor)
        hsubRun = fullfile(entries(i).paths.hsub, entries(i).fileStem);
        resFile = fullfile(hsubRun, 'TimeSeries', 'Residual', ...
            sprintf('Residual_Selective_%s_%s.xlsx', armKey, entries(i).fileStem));
        if ~isfile(resFile)
            warning('Residual file missing for %s (re-run harmonic subtraction).', entries(i).fileStem);
            continue;
        end
        if ~isfile(entries(i).paths.raw)
            warning('RAW file missing for %s; skipped.', entries(i).fileStem);
            continue;
        end

        [rawTbl, rawNames] = lb_read_excel(entries(i).paths.raw);
        resTbl = readtable(resFile);
        anchorMap = lb_anchor_map(entries(i).hsub.anchor);
        FB = lb_make_filterbank(height(resTbl), cfg.samplingMinutes, cfg);

        groups = entries(i).groups;
        for g = 1:numel(groups)
            colNames = {};
            if isfield(groups(g), 'colNames') && ~isempty(groups(g).colNames)
                colNames = groups(g).colNames;
                if isstring(colNames), colNames = cellstr(colNames); end
            elseif isfield(groups(g), 'colIdx') && ~isempty(groups(g).colIdx)
                colNames = entries(i).meta.varNames(groups(g).colIdx);
            end
            grpName = groups(g).name;

            for m = 1:numel(colNames)
                colName = char(string(colNames{m}));
                % Skip mice that failed the circadian anchor gate
                anchor = lb_anchor_lookup(anchorMap, colName);
                if ~anchor.ok, continue; end

                resVar = matlab.lang.makeValidName(['Residual_Selective_' colName]);
                if ~ismember(resVar, resTbl.Properties.VariableNames)
                    continue;
                end

                % UR measures from residual CWT
                residual = lb_prepare_signal(resTbl.(resVar));
                [wt, periods_hours] = lb_cwt(residual, FB);
                powerSpec = mean(abs(wt).^2, 2);
                bpStruct = lb_band_power(powerSpec, periods_hours, cfg);

                % CR from RAW (not residual) — circadian power on the unsubtracted signal
                crPower = NaN;
                rawColIdx = find(strcmp(rawNames, colName), 1);
                if ~isempty(rawColIdx)
                    rawSig = lb_prepare_signal(rawTbl{:, rawColIdx});
                    [wtRaw, perRaw] = lb_cwt(rawSig, FB);
                    rawBands = lb_band_power(mean(abs(wtRaw).^2, 2), perRaw, cfg);
                    crPower = rawBands.CR_20_28;
                end

                signalLabel = sprintf('%s_%s', grpName, colName);
                bpRow = struct2table(bpStruct);
                bpRow.SignalID = string(signalLabel);
                bpRow.Group = string(grpName);
                bpRow.MouseColumn = string(colName);
                bpRow.LightDuration_h = entries(i).lightHours;
                bpRow.CR_20_28 = crPower;   % overwrite residual CR with RAW CR
                bpRow.AnchorOK = true;
                bpRow.P0_h = anchor.P0_h;
                bpRow.HSubValidated = true;
                bpRow.AnalysisNote = "UR bands from selective Min360 residual; CR from RAW";
                bandRows{end+1} = bpRow; 

                % UR period peaks with harmonic-of-P0 filter
                [peaks, peakLog] = lb_validated_ur_peaks(powerSpec, periods_hours, anchor.P0_h, cfg);
                if height(peaks) > 0
                    peaks.SignalID = repmat(string(signalLabel), height(peaks), 1);
                    peaks.Group = repmat(string(grpName), height(peaks), 1);
                    peaks.MouseColumn = repmat(string(colName), height(peaks), 1);
                    peaks.LightDuration_h = repmat(entries(i).lightHours, height(peaks), 1);
                    peaks.P0_h = repmat(anchor.P0_h, height(peaks), 1);
                    peaks.HSubValidated = true(height(peaks), 1);
                    periodRows{end+1} = peaks; 
                end
                if height(peakLog) > 0
                    peakLog.SignalID = repmat(string(signalLabel), height(peakLog), 1);
                    peakLog.Group = repmat(string(grpName), height(peakLog), 1);
                    peakLog.MouseColumn = repmat(string(colName), height(peakLog), 1);
                    peakLog.LightDuration_h = repmat(entries(i).lightHours, height(peakLog), 1);
                    peakLog.P0_h = repmat(anchor.P0_h, height(peakLog), 1);
                    logRows{end+1} = peakLog; 
                end
            end
        end
    end

    bandTable = lb_vertcat_tables(bandRows);
    periodTable = lb_vertcat_tables(periodRows);
    validationLog = lb_vertcat_tables(logRows);
    if isempty(bandTable), bandTable = table(); end
    if isempty(periodTable), periodTable = table(); end
    if isempty(validationLog), validationLog = table(); end
end

function ratioTable = lb_coexpression(bandTable)
% Build CR/UR ratios used in co-expression figures.
% Main output used for plots: CR_to_UR_total = CR_20_28 / sum(UR bands).
    if isempty(bandTable) || height(bandTable) == 0
        ratioTable = table();
        return;
    end
    urBands = {'UR_1_3', 'UR_3_6', 'UR_6_12', 'UR_12_18'};
    keepCols = intersect({'SignalID', 'Group', 'LightDuration_h', 'CR_20_28'}, ...
        bandTable.Properties.VariableNames, 'stable');
    ratioTable = bandTable(:, keepCols);
    cr = bandTable.CR_20_28;
    urTot = zeros(height(bandTable), 1);
    for i = 1:numel(urBands)
        if ~ismember(urBands{i}, bandTable.Properties.VariableNames)
            continue;
        end
        ur = bandTable.(urBands{i});
        urTot = urTot + ur;
        ratioTable.(['CR_to_' urBands{i}]) = cr ./ (ur + eps);   % eps avoids /0
    end
    ratioTable.UR_total = urTot;
    ratioTable.CR_to_UR_total = cr ./ (urTot + eps);
end

%% --------------------- Wavelet maths ------------------------------------

function FB = lb_make_filterbank(nSamples, TsMinutes, cfg)
% Build the Morlet wavelet filter bank once per recording.
% Think of it as the set of "period lenses" the CWT will use (60–1680 min).
    FB = cwtfilterbank('SignalLength', nSamples, ...
        'SamplingPeriod', minutes(TsMinutes), ...
        'PeriodLimits', minutes(cfg.periodLimitsMinutes), ...
        'Wavelet', cfg.waveletName);
end

function [wt, periods_hours] = lb_cwt(signal, FB)
% Run the continuous wavelet transform.
%   wt            = complex coefficients (rows = periods, columns = time)
%   periods_hours = period (in hours) for each row of wt
% Non-finite samples are set to 0 so cwt does not fail.
    signal = signal(:);
    signal(~isfinite(signal)) = 0;
    [wt, periods] = cwt(signal, 'FilterBank', FB);
    periods_hours = hours(periods);   % convert MATLAB duration → numeric hours
end

function bandPower = lb_band_power(powerSpec, periods_hours, cfg)
% Average power inside each named band (CR 20–28, UR 1–3, …).
% powerSpec is the time-averaged |wt|² spectrum (one value per period).
    names = cfg.bands.coexpression;
    bandPower = struct();
    for i = 1:numel(names)
        bn = names{i};
        lim = cfg.bands.(bn);   % [low high] hours for this band
        mask = periods_hours >= lim(1) & periods_hours <= lim(2);
        if any(mask)
            bandPower.(bn) = mean(powerSpec(mask), 'omitnan');
        else
            bandPower.(bn) = NaN;
        end
    end
end

function peaks = lb_find_peaks(powerSpec, periods_hours, topN)
% Find the strongest peaks on log10(mean power).
% PeakRank 1 = strongest peak. Returns an empty table if none found.
    logP = log10(powerSpec + eps);
    [pks, locs, w, prom] = findpeaks(logP, 'MinPeakProminence', 0);
    if isempty(pks)
        peaks = table();
        return;
    end
    [~, ord] = sort(pks, 'descend');          % strongest first
    ord = ord(1:min(topN, numel(ord)));       % keep at most topN
    peaks = table((1:numel(ord))', periods_hours(locs(ord)), pks(ord), prom(ord), w(ord), ...
        'VariableNames', {'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', 'PeakProminence', 'PeakWidth'});
end

function [peaks, logTable] = lb_validated_ur_peaks(powerSpec, periods_hours, P0_h, cfg)
% Ultradian peaks (1–18 h) after dropping those that look like P0/k harmonics.
%   peaks    = kept peaks (ranked)
%   logTable = all candidates, with IsHarmonicOfP0 flag (for the validation log)
    urLim = cfg.hsubValidation.urPeriodRangeHours;
    urMask = periods_hours >= urLim(1) & periods_hours <= urLim(2);
    logP = log10(powerSpec + eps);
    logP_ur = logP;
    logP_ur(~urMask) = -inf;   % ignore peaks outside the UR window

    [pks, locs, w, prom] = findpeaks(logP_ur, 'MinPeakProminence', 0);
    peaks = table();
    logTable = table();
    if isempty(pks), return; end

    candPeriod = periods_hours(locs);
    isHarm = false(size(candPeriod));
    for i = 1:numel(candPeriod)
        isHarm(i) = lb_is_harmonic(candPeriod(i), P0_h, cfg);
    end
    % Record every candidate (kept and rejected) for transparency
    logTable = table(candPeriod, pks, prom, w, isHarm, ...
        'VariableNames', {'CandidatePeriod_h', 'PeakValue_log10', 'PeakProminence', ...
        'PeakWidth', 'IsHarmonicOfP0'});

    keep = ~isHarm;   % drop circadian harmonics of P0
    if ~any(keep), return; end
    pks = pks(keep); locs = locs(keep); prom = prom(keep); w = w(keep);
    candPeriod = candPeriod(keep);
    [~, ord] = sort(pks, 'descend');
    ord = ord(1:min(cfg.hsubValidation.topNPeaks, numel(ord)));
    peaks = table((1:numel(ord))', candPeriod(ord), pks(ord), prom(ord), w(ord), ...
        'VariableNames', {'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', ...
        'PeakProminence', 'PeakWidth'});
end

function isHarm = lb_is_harmonic(period_h, P0_h, cfg)
% True if period_h is close to P0/k for some small integer k
% (within ±tol, default 5%). Those peaks are treated as circadian harmonics.
    isHarm = false;
    if ~isfinite(period_h) || ~isfinite(P0_h) || P0_h <= 0, return; end
    tol = cfg.hsubValidation.harmonicToleranceFrac;
    urLim = cfg.hsubValidation.urPeriodRangeHours;
    for k = 1:cfg.hsubValidation.maxHarmonicOrder
        target = P0_h / k;   % expected period of the k-th harmonic of P0
        if target < urLim(1) || target > urLim(2), continue; end
        if abs(period_h - target) / target <= tol
            isHarm = true;
            return;
        end
    end
end

%% --------------------- Plotting (simple, readable) ----------------------

function theme = lb_theme(plotMode, cfg)
% Pick export format/DPI: development = quick PNG; publication = high-res JPEG.
    if strcmpi(plotMode, 'publication')
        theme = cfg.plot.publication;
    else
        theme = cfg.plot.development;
    end
    theme.colormap = cfg.colormap;
end

function lb_plot_scalogram(wt, periods_hours, time_day, condChangeIdx, titleStr, outPath, theme, cfg)
% Draw and save a jet time–period image (|wt|).
% White dotted lines mark light-duration changes if any were found.
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 640]);
    pcolor(time_day, periods_hours, abs(wt));   % colour = wavelet magnitude
    shading interp;
    colormap(theme.colormap);
    colorbar;
    ax = gca;
    set(ax, 'YTick', cfg.scalogramYTicksHours, 'FontName', 'Arial', 'FontSize', 11);
    hold on;
    for k = 1:numel(condChangeIdx)
        row = condChangeIdx(k) + 1;
        if row >= 1 && row <= numel(time_day)
            plot([time_day(row) time_day(row)], [min(periods_hours) max(periods_hours)], ...
                'w:', 'LineWidth', 1.5);
        end
    end
    hold off;
    xlabel('Time (days)');
    ylabel('Period (hr)');
    title(titleStr, 'Interpreter', 'none');
    lb_export_fig(fig, outPath, theme);
    close(fig);
end

function figurePaths = lb_hsub_wavelet_scalograms(hsubRoot, fileStem, groups, varNames, ...
        hsubVal, time_day, condChangeIdx, LightDur_h, FB, theme, cfg, scaloRoot)
% CWT on selective Min360 Removed and Residual from LB_01.
% Saves two-panel (Removed | Residual) group averages and individuals.
% Group averages use AnchorOK mice only.
    figurePaths = {};
    if ~isfield(hsubVal, 'available') || ~hsubVal.available
        fprintf('  Skipping HSub Removed/Residual scalograms (no HSub summary).\n');
        return;
    end

    armKey = cfg.hsubArm;   % 'Min360'
    armLabel = 'SEL_P360';
    resFile = fullfile(hsubRoot, fileStem, 'TimeSeries', 'Residual', ...
        sprintf('Residual_Selective_%s_%s.xlsx', armKey, fileStem));
    remFile = fullfile(hsubRoot, fileStem, 'TimeSeries', 'Removed', ...
        sprintf('Removed_Selective_%s_%s.xlsx', armKey, fileStem));
    if ~isfile(resFile) || ~isfile(remFile)
        warning('HSub residual/removed files missing for %s; skipped HSub scalograms.', fileStem);
        return;
    end

    resTbl = readtable(resFile);
    remTbl = readtable(remFile);
    anchorMap = lb_anchor_map(hsubVal.anchor);

    hsubScaloRoot = fullfile(scaloRoot, ['HSub_' armLabel]);
    lb_ensure_dir(hsubScaloRoot);
    fprintf('HSub Removed|Residual scalograms (%s):\n', armLabel);

    for g = 1:numel(groups)
        if isempty(groups(g).colIdx), continue; end
        grpName = groups(g).name;
        colNames = groups(g).colNames;
        if isstring(colNames), colNames = cellstr(colNames); end
        if isempty(colNames)
            colNames = varNames(groups(g).colIdx);
        end

        % Collect AnchorOK mice that exist in both residual and removed tables
        okNames = {};
        remCols = {};
        resCols = {};
        for m = 1:numel(colNames)
            colName = char(string(colNames{m}));
            anchor = lb_anchor_lookup(anchorMap, colName);
            if ~anchor.ok, continue; end
            remVar = matlab.lang.makeValidName(['Removed_Selective_' colName]);
            resVar = matlab.lang.makeValidName(['Residual_Selective_' colName]);
            if ~ismember(remVar, remTbl.Properties.VariableNames), continue; end
            if ~ismember(resVar, resTbl.Properties.VariableNames), continue; end
            okNames{end+1} = colName; %#ok<AGROW>
            remCols{end+1} = remVar; %#ok<AGROW>
            resCols{end+1} = resVar; %#ok<AGROW>
        end
        nOk = numel(okNames);
        if nOk == 0
            warning('Group "%s": no AnchorOK mice with HSub series; skipped.', grpName);
            continue;
        end

        grpFolder = fullfile(hsubScaloRoot, ['Group_' lb_sanitise(grpName)]);
        lb_ensure_dir(grpFolder);

        % --- Group average (mean across AnchorOK mice) ---
        remMat = zeros(height(remTbl), nOk);
        resMat = zeros(height(resTbl), nOk);
        for m = 1:nOk
            remMat(:, m) = lb_prepare_signal(remTbl.(remCols{m}));
            resMat(:, m) = lb_prepare_signal(resTbl.(resCols{m}));
        end
        [wtRem, periods_hours] = lb_cwt(mean(remMat, 2), FB);
        [wtRes, ~] = lb_cwt(mean(resMat, 2), FB);

        outAvg = fullfile(grpFolder, sprintf('HSub_Removed_Residual_%s_%s_L%g.%s', ...
            armLabel, lb_sanitise(grpName), LightDur_h, theme.format));
        titleAvg = sprintf('HSub %s | Average | %s | n=%d AnchorOK | L=%g h', ...
            armLabel, grpName, nOk, LightDur_h);
        lb_plot_two_panel_scalogram(wtRem, wtRes, periods_hours, time_day, ...
            condChangeIdx, titleAvg, outAvg, theme, cfg);
        figurePaths{end+1} = outAvg; %#ok<AGROW>
        fprintf('  Group average "%s" (n=%d AnchorOK) → %s\n', grpName, nOk, outAvg);

        % --- Individuals (AnchorOK only) ---
        indivFolder = fullfile(grpFolder, 'Individuals');
        lb_ensure_dir(indivFolder);
        for m = 1:nOk
            [wtRemI, ~] = lb_cwt(remMat(:, m), FB);
            [wtResI, ~] = lb_cwt(resMat(:, m), FB);
            outInd = fullfile(indivFolder, sprintf('HSub_Removed_Residual_%s_%s_%s.%s', ...
                armLabel, lb_sanitise(grpName), lb_sanitise(okNames{m}), theme.format));
            titleInd = sprintf('HSub %s | %s | %s | L=%g h', ...
                armLabel, grpName, okNames{m}, LightDur_h);
            lb_plot_two_panel_scalogram(wtRemI, wtResI, periods_hours, time_day, ...
                condChangeIdx, titleInd, outInd, theme, cfg);
            figurePaths{end+1} = outInd; %#ok<AGROW>
        end
    end
end

function lb_plot_two_panel_scalogram(wtRem, wtRes, periods_hours, time_day, ...
        condChangeIdx, titleStr, outPath, theme, cfg)
% Side-by-side Removed (left) and Residual (right) jet scalograms.
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1600 640]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    yLo = min(periods_hours); yHi = max(periods_hours);

    ax1 = nexttile(tl, 1);
    lb_draw_scalogram_axes(ax1, wtRem, periods_hours, time_day, condChangeIdx, ...
        yLo, yHi, 'Removed', theme, cfg);
    ax2 = nexttile(tl, 2);
    lb_draw_scalogram_axes(ax2, wtRes, periods_hours, time_day, condChangeIdx, ...
        yLo, yHi, 'Residual', theme, cfg);

    title(tl, titleStr, 'Interpreter', 'none', 'FontName', 'Arial');
    lb_export_fig(fig, outPath, theme);
    close(fig);
end

function lb_draw_scalogram_axes(ax, wt, periods_hours, time_day, condChangeIdx, ...
        yLo, yHi, panelLabel, theme, cfg)
% One panel of a jet scalogram (used by single- and two-panel exports).
    pcolor(ax, time_day, periods_hours, abs(wt));
    shading(ax, 'interp');
    colormap(ax, theme.colormap);
    colorbar(ax);
    set(ax, 'YTick', cfg.scalogramYTicksHours, 'FontName', 'Arial', 'FontSize', 11);
    ylim(ax, [yLo yHi]);
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
    title(ax, panelLabel, 'Interpreter', 'none');
end

function lb_plot_coexpression(ratioTable, outPath, titleStr, ylabelStr)
% Line plot: mean CR/UR_total vs light duration, one line per condition group.
    if isempty(ratioTable) || height(ratioTable) == 0
        warning('No co-expression data to plot: %s', outPath);
        return;
    end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 500]);
    groups = unique(ratioTable.Group);
    hold on;
    cols = lines(max(numel(groups), 1));
    for g = 1:numel(groups)
        sub = ratioTable(ratioTable.Group == groups(g), :);
        x = unique(sub.LightDuration_h);
        mu = zeros(size(x));
        for k = 1:numel(x)
            % Average across mice at this light duration
            mu(k) = mean(sub.CR_to_UR_total(sub.LightDuration_h == x(k)), 'omitnan');
        end
        plot(x, mu, '-o', 'Color', cols(g,:), 'LineWidth', 1.5, ...
            'DisplayName', char(groups(g)));
    end
    hold off;
    xlabel('Light duration (h)');
    ylabel(ylabelStr, 'Interpreter', 'tex');
    title(titleStr);
    legend('Location', 'best');
    set(gca, 'FontName', 'Arial', 'FontSize', 11);
    [~, ~, ext] = fileparts(outPath);
    theme = struct('format', erase(ext, '.'), 'dpi', 96);
    if strcmpi(theme.format, 'jpeg'), theme.dpi = 600; end
    lb_export_fig(fig, outPath, theme);
    close(fig);
end

function lb_plot_period_scatter(periodTable, outPath, titleStr, ylabelStr, ylimRange)
% Scatter of rank-1 peak period vs light duration (one point per mouse/average).
    if isempty(periodTable) || height(periodTable) == 0
        warning('No period data to plot: %s', outPath);
        return;
    end
    if ~ismember('PeakRank', periodTable.Properties.VariableNames)
        return;
    end
    top1 = periodTable(periodTable.PeakRank == 1, :);   % strongest peak only
    if isempty(top1), return; end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 500]);
    groups = unique(top1.Group);
    cols = lines(max(numel(groups), 1));
    hold on;
    for g = 1:numel(groups)
        sub = top1(top1.Group == groups(g), :);
        scatter(sub.LightDuration_h, sub.PeakPeriod_hr, 36, cols(g,:), 'filled', ...
            'DisplayName', char(groups(g)));
    end
    hold off;
    xlabel('Light duration (h)');
    ylabel(ylabelStr);
    title(titleStr);
    ylim(ylimRange);
    legend('Location', 'best');
    set(gca, 'FontName', 'Arial', 'FontSize', 11);
    [~, ~, ext] = fileparts(outPath);
    theme = struct('format', erase(ext, '.'), 'dpi', 96);
    if strcmpi(theme.format, 'jpeg'), theme.dpi = 600; end
    lb_export_fig(fig, outPath, theme);
    close(fig);
end

function lb_export_fig(fig, outPath, theme)
% Save the figure to disk (PNG or JPEG). Falls back to saveas if exportgraphics fails.
    [folder, ~, ~] = fileparts(outPath);
    if ~isempty(folder), lb_ensure_dir(folder); end
    try
        exportgraphics(fig, outPath, 'Resolution', theme.dpi);
    catch
        saveas(fig, outPath);
    end
end

%% --------------------- Across-photoperiod extra figures -----------------
% Period distributions and stitched scalograms complete the cohort figure set.
% Distributions = one box+jitter panel per condition group.
% Stitched = concatenate photoperiod segments (ordered by light duration),
% then compute one CWT on the long signal (RAW or Removed|Residual).

function outPaths = lb_plot_period_distributions(periodTable, figFolder, baseName, theme, ...
        titlePrefix, ylabelStr, ylimRange)
% One box+jitter figure per condition group (rank-1 peaks only).
% Box = median/IQR at each light duration; points = individual mice (jittered).
    outPaths = {};
    if isempty(periodTable) || height(periodTable) == 0
        warning('No period data for distributions: %s', baseName);
        return;
    end
    if ~ismember('PeakRank', periodTable.Properties.VariableNames)
        return;
    end
    top1 = periodTable(periodTable.PeakRank == 1, :);   % strongest peak only
    if isempty(top1), return; end

    distFolder = fullfile(figFolder, 'Period_Distributions');
    lb_ensure_dir(distFolder);
    groups = unique(top1.Group);
    cols = lines(max(numel(groups), 1));

    for g = 1:numel(groups)
        sub = top1(top1.Group == groups(g), :);
        if isempty(sub), continue; end
        outFile = fullfile(distFolder, sprintf('%s_%s.%s', ...
            baseName, lb_sanitise(char(groups(g))), theme.format));

        photos = sort(unique(sub.LightDuration_h));
        photos = photos(isfinite(photos));
        if isempty(photos), continue; end

        % Stack (x,y) for boxchart: x = light duration, y = period (h)
        xAll = [];
        yAll = [];
        for k = 1:numel(photos)
            yk = sub.PeakPeriod_hr(sub.LightDuration_h == photos(k));
            yk = yk(isfinite(yk));
            xAll = [xAll; repmat(photos(k), numel(yk), 1)]; %#ok<AGROW>
            yAll = [yAll; yk(:)]; %#ok<AGROW>
        end
        if isempty(yAll), continue; end

        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 960 520]);
        ax = axes(fig);
        hold(ax, 'on');
        boxchart(ax, xAll, yAll, ...
            'BoxFaceColor', cols(g,:), 'BoxFaceAlpha', 1, ...
            'BoxEdgeColor', [0 0 0], 'LineWidth', 1, ...
            'MarkerStyle', 'none', 'BoxWidth', 0.55);
        rng(42);   % fixed jitter so re-runs look the same
        for k = 1:numel(photos)
            yk = sub.PeakPeriod_hr(sub.LightDuration_h == photos(k));
            yk = yk(isfinite(yk));
            n = numel(yk);
            if n == 0, continue; end
            jitter = (rand(n, 1) - 0.5) * 0.7;
            scatter(ax, photos(k) + jitter, yk, 36, cols(g,:), 'filled');
        end
        hold(ax, 'off');
        xticks(ax, photos);
        xlabel(ax, 'Light duration (h)');
        ylabel(ax, ylabelStr);
        title(ax, sprintf('%s — %s', titlePrefix, char(groups(g))), 'Interpreter', 'none');
        ylim(ax, ylimRange);
        if numel(photos) >= 2
            pad = 0.15 * (photos(end) - photos(1));
            xlim(ax, [photos(1) - pad, photos(end) + pad]);
        end
        set(ax, 'FontName', 'Arial', 'FontSize', 11);
        lb_export_fig(fig, outFile, theme);
        close(fig);
        outPaths{end+1} = outFile; %#ok<AGROW>
        fprintf('  Period distribution: %s → %s\n', char(groups(g)), outFile);
    end
end

function figurePaths = lb_plot_stitched_scalograms(entries, figFolder, cohortTag, theme, cfg)
% Group-average RAW activity stitched across photoperiods → one scalogram.
% White lines mark photoperiod joins; L-tags sit above each segment.
    figurePaths = {};
    groupNames = lb_collect_group_names(entries);
    if isempty(groupNames)
        warning('No condition groups for stitched RAW scalograms.');
        return;
    end
    stitchedFolder = fullfile(figFolder, 'Stitched_Scalograms');
    lb_ensure_dir(stitchedFolder);

    for g = 1:numel(groupNames)
        grpName = groupNames{g};
        stitch = lb_stitch_raw_group(entries, grpName, cfg);
        if isempty(stitch.signal)
            warning('Skipping stitched RAW for "%s": no valid segments.', grpName);
            continue;
        end
        FB = lb_make_filterbank(numel(stitch.signal), cfg.samplingMinutes, cfg);
        [wt, periods_hours] = lb_cwt(stitch.signal, FB);
        outFile = fullfile(stitchedFolder, sprintf('Stitched_Scalogram_Average_%s_%s.%s', ...
            lb_sanitise(grpName), cohortTag, theme.format));
        titleStr = sprintf('RAW | Stitched average | %s | %s | L%g–L%g h', ...
            grpName, cohortTag, min(stitch.photoHours), max(stitch.photoHours));
        lb_render_stitched_scalogram(wt, periods_hours, stitch, titleStr, outFile, theme, cfg);
        figurePaths{end+1} = outFile; %#ok<AGROW>
        fprintf('  Stitched RAW scalogram: "%s" → %s\n', grpName, outFile);
    end
end

function figurePaths = lb_plot_stitched_hsub_scalograms(entries, figFolder, cohortTag, theme, cfg)
% Stitched SEL_P360 Removed|Residual group averages across photoperiod.
% Reads selective Min360 Removed and Residual workbooks from harmonic subtraction.
    figurePaths = {};
    groupNames = lb_collect_group_names(entries);
    if isempty(groupNames)
        warning('No condition groups for stitched HSub scalograms.');
        return;
    end
    stitchedFolder = fullfile(figFolder, 'Stitched_Scalograms_HSub');
    lb_ensure_dir(stitchedFolder);
    armKey = cfg.hsubArm;
    armLabel = 'SEL_P360';

    for g = 1:numel(groupNames)
        grpName = groupNames{g};
        [stitchRem, stitchRes] = lb_stitch_hsub_group(entries, grpName, armKey, cfg);
        if isempty(stitchRem.signal) || isempty(stitchRes.signal)
            warning('Skipping stitched HSub for "%s": no valid segments.', grpName);
            continue;
        end
        FB = lb_make_filterbank(numel(stitchRem.signal), cfg.samplingMinutes, cfg);
        [wtRem, periods_hours] = lb_cwt(stitchRem.signal, FB);
        [wtRes, ~] = lb_cwt(stitchRes.signal, FB);
        outFile = fullfile(stitchedFolder, sprintf('Stitched_HSub_Removed_Residual_%s_%s_%s.%s', ...
            armLabel, lb_sanitise(grpName), cohortTag, theme.format));
        titleStr = sprintf('HSub %s | Stitched average | %s | %s | L%g–L%g h', ...
            armLabel, grpName, cohortTag, min(stitchRem.photoHours), max(stitchRem.photoHours));
        lb_render_stitched_hsub(wtRem, wtRes, periods_hours, stitchRem, titleStr, outFile, theme, cfg);
        figurePaths{end+1} = outFile; %#ok<AGROW>
        fprintf('  Stitched HSub scalogram: "%s" → %s\n', grpName, outFile);
    end
end

function names = lb_collect_group_names(entries)
% Unique condition-group names across the cohort (stable order of first appearance).
    names = {};
    for i = 1:numel(entries)
        if ~isfield(entries(i), 'groups') || isempty(entries(i).groups), continue; end
        for g = 1:numel(entries(i).groups)
            names{end+1} = entries(i).groups(g).name; %#ok<AGROW>
        end
    end
    names = unique(names, 'stable');
end

function stitch = lb_stitch_raw_group(entries, grpName, cfg)
% Build a long RAW group-mean signal by concatenating each photoperiod segment.
% Time axis is continuous in "days" so the CWT sees one unbroken series.
    stitch = lb_init_stitch();
    dayOffset = 0;
    for i = 1:numel(entries)
        grp = lb_find_group(entries(i).groups, grpName);
        if isempty(grp), continue; end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'raw')
            warning('paths.raw missing for %s; segment skipped.', entries(i).fileStem);
            continue;
        end
        if ~isfile(entries(i).paths.raw)
            warning('RAW file missing for %s; segment skipped.', entries(i).fileStem);
            continue;
        end
        [tbl, varNames] = lb_read_excel(entries(i).paths.raw);
        % Map Time/Light quietly (no mouse-column dialog — batch must stay unattended)
        [~, lightIdx] = lb_map_time_light_only(varNames);
        colIdx = lb_group_col_idx(grp, varNames);
        if isempty(colIdx), continue; end

        segSignal = mean(lb_stack_signals(tbl, colIdx), 2);
        n = numel(segSignal);
        segDay = ((0:n-1)' * cfg.samplingMinutes) / (60 * 24) + dayOffset;
        lightBounds = [];
        if ~isempty(lightIdx)
            lightVec = tbl{:, lightIdx};
            if ~isnumeric(lightVec), lightVec = str2double(string(lightVec)); end
            condChangeIdx = find(diff(lightVec) ~= 0);
            for k = 1:numel(condChangeIdx)
                row = condChangeIdx(k) + 1;
                if row >= 1 && row <= n
                    lightBounds(end+1) = segDay(row); %#ok<AGROW>
                end
            end
        end
        stitch = lb_append_segment(stitch, segSignal, segDay, dayOffset, ...
            entries(i).lightHours, lightBounds, cfg);
        dayOffset = stitch.nextDayOffset;
    end
end

function [stitchRem, stitchRes] = lb_stitch_hsub_group(entries, grpName, armKey, cfg)
% Same stitching idea as RAW, but using selective Min360 Removed and Residual.
    stitchRem = lb_init_stitch();
    stitchRes = lb_init_stitch();
    dayOffset = 0;
    for i = 1:numel(entries)
        grp = lb_find_group(entries(i).groups, grpName);
        if isempty(grp), continue; end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'hsub')
            warning('paths.hsub missing for %s; segment skipped.', entries(i).fileStem);
            continue;
        end
        hsubRun = fullfile(entries(i).paths.hsub, entries(i).fileStem);
        [avgRem, avgRes, segDay, lightBounds] = lb_load_hsub_group_avg( ...
            hsubRun, entries(i).fileStem, grp, armKey, cfg);
        if isempty(avgRem), continue; end
        segDay = segDay + dayOffset;
        stitchRem = lb_append_segment(stitchRem, avgRem, segDay, dayOffset, ...
            entries(i).lightHours, lightBounds, cfg);
        stitchRes = lb_append_segment(stitchRes, avgRes, segDay, dayOffset, ...
            entries(i).lightHours, lightBounds, cfg);
        dayOffset = stitchRem.nextDayOffset;
    end
end

function [avgRem, avgRes, time_day, lightBounds] = lb_load_hsub_group_avg( ...
        hsubRunFolder, fileStem, grp, armKey, cfg)
% Mean Removed and Residual series for one group within one photoperiod file.
    avgRem = []; avgRes = []; time_day = []; lightBounds = [];
    remFile = fullfile(hsubRunFolder, 'TimeSeries', 'Removed', ...
        sprintf('Removed_Selective_%s_%s.xlsx', armKey, fileStem));
    resFile = fullfile(hsubRunFolder, 'TimeSeries', 'Residual', ...
        sprintf('Residual_Selective_%s_%s.xlsx', armKey, fileStem));
    if ~isfile(remFile) || ~isfile(resFile)
        warning('HSub timeseries missing for %s.', fileStem);
        return;
    end
    remTbl = readtable(remFile);
    resTbl = readtable(resFile);
    colNames = {};
    if isfield(grp, 'colNames') && ~isempty(grp.colNames)
        colNames = grp.colNames;
        if isstring(colNames), colNames = cellstr(colNames); end
    end
    remCols = [];
    resCols = [];
    for c = 1:numel(colNames)
        remVar = matlab.lang.makeValidName(['Removed_Selective_' colNames{c}]);
        resVar = matlab.lang.makeValidName(['Residual_Selective_' colNames{c}]);
        if ismember(remVar, remTbl.Properties.VariableNames)
            remCols(:, end+1) = remTbl.(remVar); %#ok<AGROW>
        end
        if ismember(resVar, resTbl.Properties.VariableNames)
            resCols(:, end+1) = resTbl.(resVar); %#ok<AGROW>
        end
    end
    if isempty(remCols) || isempty(resCols), return; end
    avgRem = mean(remCols, 2, 'omitnan');
    avgRes = mean(resCols, 2, 'omitnan');
    n = numel(avgRem);
    time_day = ((0:n-1)' * cfg.samplingMinutes) / (60 * 24);
    lightVar = remTbl.Properties.VariableNames{end};
    if contains(lower(lightVar), 'light')
        lightVec = remTbl.(lightVar);
        condChangeIdx = find(diff(lightVec) ~= 0);
        for k = 1:numel(condChangeIdx)
            row = condChangeIdx(k) + 1;
            if row >= 1 && row <= n
                lightBounds(end+1) = time_day(row); %#ok<AGROW>
            end
        end
    end
end

function stitch = lb_init_stitch()
% Empty stitch struct: signal, continuous time axis, and boundary markers.
    stitch = struct('signal', [], 'time_day', [], 'photoBounds', [], ...
        'lightBounds', [], 'photoHours', [], 'segmentMid', [], 'nextDayOffset', 0);
end

function stitch = lb_append_segment(stitch, segSignal, segDay, dayOffset, ...
        lightHours, lightBounds, cfg)
% Append one photoperiod segment onto the growing stitch.
% photoBounds mark joins between photoperiods; lightBounds mark within-file
% light-duration changes (usually none when one photoperiod per file).
    TsDay = cfg.samplingMinutes / (60 * 24);
    segStart = dayOffset;
    if ~isempty(stitch.signal)
        stitch.photoBounds(end+1) = dayOffset; %#ok<AGROW>
    end
    if ~isempty(lightBounds)
        stitch.lightBounds = [stitch.lightBounds, lightBounds(:)']; %#ok<AGROW>
    end
    stitch.signal = [stitch.signal; segSignal(:)]; %#ok<AGROW>
    stitch.time_day = [stitch.time_day; segDay(:)]; %#ok<AGROW>
    stitch.photoHours(end+1) = lightHours; %#ok<AGROW>
    stitch.segmentMid(end+1) = (segStart + segDay(end)) / 2; %#ok<AGROW>
    stitch.nextDayOffset = segDay(end) + TsDay;
end

function grp = lb_find_group(groups, grpName)
% Return the group struct with matching name, or [] if not found.
    grp = [];
    for k = 1:numel(groups)
        if strcmp(groups(k).name, grpName)
            grp = groups(k);
            return;
        end
    end
end

function colIdx = lb_group_col_idx(grp, varNames)
% Resolve group membership to column indices via stored mouse names.
    colIdx = [];
    if isfield(grp, 'colNames') && ~isempty(grp.colNames)
        colNames = grp.colNames;
        if isstring(colNames), colNames = cellstr(colNames); end
        for c = 1:numel(colNames)
            ix = find(strcmp(varNames, colNames{c}), 1);
            if ~isempty(ix), colIdx(end+1) = ix; end %#ok<AGROW>
        end
    elseif isfield(grp, 'colIdx')
        colIdx = grp.colIdx;
    end
end

function [timeIdx, lightIdx] = lb_map_time_light_only(varNames)
% Map Time / Light columns without prompting for mouse columns (batch-safe).
    timeIdx = []; lightIdx = [];
    for i = 1:numel(varNames)
        vn = lower(strtrim(varNames{i}));
        if contains(vn, 'time') && isempty(timeIdx), timeIdx = i; end
        if contains(vn, 'light') && contains(vn, 'duration'), lightIdx = i; end
    end
end

function lb_render_stitched_scalogram(wt, periods_hours, stitch, titleStr, outPath, theme, cfg)
% Single-panel jet scalogram of a stitched RAW series.
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1280 640]);
    pcolor(stitch.time_day, periods_hours, abs(wt));
    shading interp;
    colormap(theme.colormap);
    colorbar;
    ax = gca;
    set(ax, 'YTick', cfg.scalogramYTicksHours, 'FontName', 'Arial', 'FontSize', 11);
    yLo = min(periods_hours); yHi = max(periods_hours);
    hold on;
    for k = 1:numel(stitch.photoBounds)
        plot([stitch.photoBounds(k) stitch.photoBounds(k)], [yLo yHi], 'w:', 'LineWidth', 2);
    end
    for k = 1:numel(stitch.lightBounds)
        plot([stitch.lightBounds(k) stitch.lightBounds(k)], [yLo yHi], 'w:', 'LineWidth', 1.5);
    end
    for k = 1:numel(stitch.segmentMid)
        text(stitch.segmentMid(k), yHi * 0.96, sprintf('L%g', stitch.photoHours(k)), ...
            'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 9);
    end
    hold off;
    xlabel('Time (days)');
    ylabel('Period (hr)');
    title(titleStr, 'Interpreter', 'none');
    lb_export_fig(fig, outPath, theme);
    close(fig);
end

function lb_render_stitched_hsub(wtRem, wtRes, periods_hours, stitch, titleStr, outPath, theme, cfg)
% Two-panel stitched Removed|Residual scalogram (same time axis on both).
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1680 640]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    yLo = min(periods_hours); yHi = max(periods_hours);
    lb_draw_stitched_panel(nexttile(tl, 1), wtRem, periods_hours, stitch, yLo, yHi, 'Removed', theme, cfg);
    lb_draw_stitched_panel(nexttile(tl, 2), wtRes, periods_hours, stitch, yLo, yHi, 'Residual', theme, cfg);
    sgtitle(fig, titleStr, 'Interpreter', 'none', 'FontName', 'Arial');
    lb_export_fig(fig, outPath, theme);
    close(fig);
end

function lb_draw_stitched_panel(ax, wt, periods_hours, stitch, yLo, yHi, panelLabel, theme, cfg)
% One panel of a stitched two-panel figure (photoperiod joins + L-tags).
    axes(ax); %#ok<LAXES>
    pcolor(stitch.time_day, periods_hours, abs(wt));
    shading interp;
    colormap(ax, theme.colormap);
    colorbar;
    set(ax, 'YTick', cfg.scalogramYTicksHours, 'FontName', 'Arial', 'FontSize', 11);
    hold(ax, 'on');
    for k = 1:numel(stitch.photoBounds)
        x = stitch.photoBounds(k);
        plot(ax, [x x], [yLo yHi], 'w:', 'LineWidth', 2);
    end
    for k = 1:numel(stitch.lightBounds)
        x = stitch.lightBounds(k);
        plot(ax, [x x], [yLo yHi], 'w:', 'LineWidth', 1.5);
    end
    for k = 1:numel(stitch.segmentMid)
        text(ax, stitch.segmentMid(k), yHi * 0.96, sprintf('L%g', stitch.photoHours(k)), ...
            'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 9);
    end
    hold(ax, 'off');
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Period (hr)');
    title(ax, panelLabel);
end

function lb_update_summaries_index(summaryDir, fileStem, lightHours, summaryPath)
% Append/replace one row in Summaries_Index.xlsx so the cohort list stays current.
    idxPath = fullfile(summaryDir, 'Summaries_Index.xlsx');
    row = table(string(fileStem), lightHours, string(summaryPath), datetime('now'), ...
        'VariableNames', {'FileStem', 'LightDuration_h', 'SummaryPath', 'Created'});
    if isfile(idxPath)
        try
            old = readtable(idxPath);
            keep = old.FileStem ~= string(fileStem);
            old = old(keep, :);
            common = intersect(old.Properties.VariableNames, row.Properties.VariableNames, 'stable');
            row = [old(:, common); row(:, common)];
        catch
            % overwrite on read failure
        end
    end
    writetable(row, idxPath);
end

%% --------------------- Groups / I/O / misc ------------------------------

function groups = lb_group_dialog(varNames, mouseIdx)
% Interactive condition groups (same flow as the main analysis codes):
%   1. Type how many groups
%   2. For each group: type the name, then tick which mouse columns belong to it
%
% Returns a struct array with fields: name, colIdx, colNames.
    mouseNames = varNames(mouseIdx);

    % --- Step 1: how many groups? ---
    ansCell = inputdlg({'Enter number of condition groups:'}, ...
        'Condition Groups', [1 50], {'2'});
    if isempty(ansCell)
        error('Group assignment cancelled.');
    end
    nG = str2double(strtrim(ansCell{1}));
    if isnan(nG) || nG < 1 || mod(nG, 1) ~= 0
        error('Invalid number of groups (need a positive whole number).');
    end

    % --- Step 2: name + column picker for each group ---
    groups = struct('name', {}, 'colIdx', {}, 'colNames', {});
    for g = 1:nG
        grpNameCell = inputdlg( ...
            {sprintf('Enter a name for condition group %d:', g)}, ...
            'Condition Group Name', [1 50]);
        if isempty(grpNameCell)
            error('No name entered for group %d.', g);
        end
        grpName = strtrim(grpNameCell{1});
        if isempty(grpName)
            error('Empty name for group %d.', g);
        end

        [sel, ok] = listdlg( ...
            'PromptString', sprintf('Select activity columns for group "%s":', grpName), ...
            'SelectionMode', 'multiple', ...
            'ListString', mouseNames, ...
            'ListSize', [300 400]);
        if ~ok || isempty(sel)
            error('No columns selected for group "%s".', grpName);
        end

        groups(g).name = grpName;
        groups(g).colIdx = mouseIdx(sel);      % indices into the full table
        groups(g).colNames = mouseNames(sel);  % stored for reuse on later files
        fprintf('  Group "%s": %d mice\n', grpName, numel(groups(g).colIdx));
    end

    fprintf('Assigned %d condition group(s).\n', nG);
end

function groups = lb_resolve_groups(varNames, template)
% Re-apply a previous group definition to a new file by matching mouse names.
    groups = template;
    for g = 1:numel(groups)
        if isfield(template(g), 'colNames') && ~isempty(template(g).colNames)
            names = template(g).colNames;
            if isstring(names), names = cellstr(names); end
            idx = zeros(1, 0);
            for k = 1:numel(names)
                hit = find(strcmp(varNames, names{k}), 1);
                if ~isempty(hit), idx(end+1) = hit; end %#ok<AGROW>
            end
            groups(g).colIdx = idx;
            groups(g).colNames = varNames(idx);
        end
    end
end

function M = lb_stack_signals(tbl, colIdx)
% Put several mouse columns side-by-side as a matrix (rows = time, cols = mice).
    M = zeros(height(tbl), numel(colIdx));
    for i = 1:numel(colIdx)
        M(:, i) = lb_prepare_signal(tbl{:, colIdx(i)});
    end
end

function signal = lb_prepare_signal(x)
% Make a clean column vector for CWT: numeric, non-finite → 0.
    signal = x(:);
    if ~isnumeric(signal), signal = str2double(string(signal)); end
    signal(~isfinite(signal)) = 0;
end

function map = lb_anchor_map(anchorTable)
% Turn the Anchor_Report table into a lookup: mouse name → (ok, P0_h).
    map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if isempty(anchorTable) || ~istable(anchorTable), return; end
    colVar = lb_pick_var(anchorTable, {'Column', 'colName', 'Mouse'});
    okVar = lb_pick_var(anchorTable, {'AnchorOK', 'anchorOK'});
    p0Var = lb_pick_var(anchorTable, {'Period_hours', 'P0_h', 'Period_h'});
    for r = 1:height(anchorTable)
        key = char(string(anchorTable.(colVar)(r)));
        okVal = anchorTable.(okVar)(r);
        % AnchorOK may be logical, numeric, or text depending on how Excel was saved
        if islogical(okVal), ok = okVal;
        elseif isnumeric(okVal), ok = okVal ~= 0;
        else, ok = any(strcmpi(string(okVal), ["true", "1", "yes"]));
        end
        p0 = anchorTable.(p0Var)(r);
        if ~isnumeric(p0), p0 = str2double(string(p0)); end
        map(key) = struct('ok', ok, 'P0_h', double(p0));
    end
end

function anchor = lb_anchor_lookup(map, colName)
% Look up one mouse in the anchor map. ok is true only if AnchorOK and P0 is valid.
    anchor = struct('ok', false, 'P0_h', NaN);
    if ~isKey(map, colName), return; end
    anchor = map(colName);
    anchor.ok = anchor.ok && isfinite(anchor.P0_h) && anchor.P0_h > 0;
end

function v = lb_pick_var(tbl, candidates)
% Return the first column name from candidates that exists in tbl.
    for i = 1:numel(candidates)
        if ismember(candidates{i}, tbl.Properties.VariableNames)
            v = candidates{i}; return;
        end
    end
    error('Expected one of: %s', strjoin(candidates, ', '));
end

function out = lb_vertcat_tables(cells)
% Stack a cell array of tables into one table (skipping empties).
    cells = cells(~cellfun(@isempty, cells));
    if isempty(cells), out = table(); return; end
    out = cells{1};
    for i = 2:numel(cells)
        out = lb_vertcat_two(out, cells{i});
    end
end

function out = lb_vertcat_two(A, B)
% Stack two tables after aligning their column names.
% (Different mice/groups can produce tables with slightly different columns.)
    if isempty(A), out = B; return; end
    if isempty(B), out = A; return; end
    allNames = unique([A.Properties.VariableNames, B.Properties.VariableNames], 'stable');
    for i = 1:numel(allNames)
        nm = allNames{i};
        if ~ismember(nm, A.Properties.VariableNames)
            A.(nm) = lb_missing_like(B.(nm), height(A));
        end
        if ~ismember(nm, B.Properties.VariableNames)
            B.(nm) = lb_missing_like(A.(nm), height(B));
        end
    end
    A = A(:, allNames);
    B = B(:, allNames);
    out = [A; B];   % vertical concatenate
end

function col = lb_missing_like(template, n)
% Make an n-row "empty" column with the same type as template (NaN / "" / false).
    if isstring(template)
        col = strings(n, 1);
    elseif iscell(template)
        col = cell(n, 1);
    elseif islogical(template)
        col = false(n, 1);
    else
        col = NaN(n, 1);
    end
end

function [fileList, inPath] = lb_select_input_files(runMode)
% Pop-up dialog(s) to choose one Excel/CSV file, or several from a folder.
    fileList = {};
    inPath = pwd;
    if strcmpi(runMode, 'Single file')
        [f, p] = uigetfile({'*.xlsx;*.xls;*.csv', 'Excel / CSV'}, 'Select activity Excel file');
        if isequal(f, 0), return; end   % Cancel
        inPath = p;
        fileList = {fullfile(p, f)};
    else
        inPath = uigetdir(pwd, 'Select folder containing activity Excel files');
        if isequal(inPath, 0), inPath = pwd; return; end
        d = [dir(fullfile(inPath, '*.xlsx')); dir(fullfile(inPath, '*.csv'))];
        d = d(~[d.isdir]);
        if isempty(d), return; end
        names = unique({d.name}, 'stable');
        [sel, ok] = listdlg('PromptString', 'Select files:', ...
            'SelectionMode', 'multiple', 'ListString', names, 'InitialValue', 1:numel(names));
        if ~ok || isempty(sel), return; end
        fileList = fullfile(inPath, names(sel));
        if ischar(fileList), fileList = {fileList}; end
    end
end

function [tbl, varNames] = lb_read_excel(inputFile)
% Read Excel/CSV into a table; keep original header spellings.
    tbl = readtable(inputFile, 'VariableNamingRule', 'preserve');
    varNames = tbl.Properties.VariableNames;
end

function [timeIdx, lightIdx, mouseIdx] = lb_map_columns(varNames)
% Find Time, Light duration, and mouse activity columns by header name.
    timeIdx = []; lightIdx = [];
    for i = 1:numel(varNames)
        vn = lower(strtrim(varNames{i}));
        if contains(vn, 'time') && isempty(timeIdx), timeIdx = i; end
        if contains(vn, 'light') && contains(vn, 'duration'), lightIdx = i; end
    end
    mouseIdx = setdiff(1:numel(varNames), [timeIdx, lightIdx]);
    keep = true(size(mouseIdx));
    for k = 1:numel(mouseIdx)
        vn = lower(varNames{mouseIdx(k)});
        if contains(vn, 'zt') || strcmp(vn, 'day') || contains(vn, 'unnamed')
            keep(k) = false;
        end
    end
    mouseIdx = mouseIdx(keep);
end

function s = lb_sanitise(strIn)
% Make a string safe for filenames (replace odd characters with _).
    s = char(regexprep(string(strIn), '[^\w\-]', '_'));
end

function lb_ensure_dir(p)
% Create folder p if missing.
    if ~exist(p, 'dir'), mkdir(p); end
end
