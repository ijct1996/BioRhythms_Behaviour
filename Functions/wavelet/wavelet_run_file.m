function summary = wavelet_run_file(rawFile, hsubFolder, outputFolder, varargin)
%WAVELET_RUN_FILE Wavelet on RAW (+ HSub validation); scalograms, peaks, handoff.

    p = inputParser;
    addParameter(p, 'Interactive', true, @islogical);
    addParameter(p, 'PlotMode', 'development', @ischar);
    addParameter(p, 'Groups', [], @(x) isempty(x) || isstruct(x));
    parse(p, varargin{:});

    theme = plot_config(p.Results.PlotMode);
    theme = plot_theme_ensure_scalogram(theme);
    cfg = core_defaults();

    [tbl, meta] = read_behav_excel(rawFile);
    if isempty(p.Results.Groups)
        groups = group_assignment_dialog(meta.varNames, meta.mouseIdx, meta, ...
            'Interactive', p.Results.Interactive);
        groups = groups_to_template(groups, meta.varNames);
    else
        groups = resolve_groups_by_names(meta.varNames, p.Results.Groups);
        fprintf('Reusing %d condition group(s) for %s.\n', numel(groups), meta.fileStem);
    end

    time_hours = tbl{:, meta.timeIdx};
    time_min = time_hours * 60;
    time_day = time_min / (60 * 24);
    TsMinutes = median(diff(time_min), 'omitnan');
    lightVec = tbl{:, meta.lightIdx};
    condChangeIdx = find(diff(lightVec) ~= 0);

    [~, fileStem, ~] = fileparts(rawFile);
    runFolder = fullfile(outputFolder, fileStem);
    scaloRoot = fullfile(runFolder, 'Scalograms');
    peakFolder = fullfile(runFolder, 'PeriodPeaks');
    handoffDir = fullfile(fileparts(outputFolder), 'Handoff');
    ensure_dir(scaloRoot);
    ensure_dir(peakFolder);
    ensure_dir(handoffDir);

    hsubSummary = fullfile(hsubFolder, fileStem, 'Reports', 'HarmonicRemoval_Summary.xlsx');
    hsubVal = struct('available', isfile(hsubSummary), 'path', hsubSummary);
    if hsubVal.available
        try
            hsubVal.anchor = readtable(hsubSummary, 'Sheet', 'Anchor_Report');
            hsubVal.recommendation = readtable(hsubSummary, 'Sheet', 'Recommendation');
        catch
            hsubVal.available = false;
        end
    end

    periodAll = {};
    bandAll = {};
    figurePaths = {};
    FB = wavelet_make_filterbank(height(tbl), TsMinutes, cfg);

    % Group-average scalogram per user-defined condition group
    fprintf('Group-average scalograms for %s:\n', fileStem);
    [avgPaths, avgPeriods, avgBands] = wavelet_plot_group_averages( ...
        tbl, meta, groups, FB, time_day, condChangeIdx, scaloRoot, theme, cfg);
    figurePaths = [figurePaths, avgPaths];
    periodAll = [periodAll, avgPeriods];
    bandAll = [bandAll, avgBands];

    % Individual scalograms within each group
    ext = theme.scalogram.format;
    for g = 1:numel(groups)
        grpName = groups(g).name;
        grpFolder = fullfile(scaloRoot, ['Group_' sanitise_filename(grpName)]);
        indivFolder = fullfile(grpFolder, 'Individuals');
        ensure_dir(indivFolder);

        for idx = 1:numel(groups(g).colIdx)
            colIdx = groups(g).colIdx(idx);
            signalID = meta.varNames{colIdx};
            signal = wavelet_prepare_signal(tbl, colIdx);

            [wt, periods_hours, ~] = wavelet_compute_cwt(signal, FB);
            powerSpec = mean(abs(wt).^2, 2);

            outScalo = fullfile(indivFolder, sprintf('Scalogram_%s_%s.%s', ...
                sanitise_filename(grpName), sanitise_filename(signalID), ext));
            wavelet_plot_scalogram(wt, periods_hours, time_day, condChangeIdx, ...
                sprintf('RAW | %s | %s', signalID, grpName), outScalo, theme);
            figurePaths{end+1} = outScalo; %#ok<AGROW>

            peaks = wavelet_find_period_peaks(powerSpec, periods_hours, cfg.wavelet.topNPeaks);
            peaks.SignalID = repmat(string(sprintf('%s_%s', grpName, signalID)), height(peaks), 1);
            peaks.Group = repmat(string(grpName), height(peaks), 1);
            peaks.LightDuration_h = repmat(meta.lightDurationHours, height(peaks), 1);
            peaks.N_Mice = ones(height(peaks), 1);
            if height(peaks) > 0
                periodAll{end+1} = peaks; %#ok<AGROW>
            end

            bp = wavelet_band_power(powerSpec, periods_hours, cfg);
            bpRow = struct2table(bp);
            bpRow.SignalID = string(sprintf('%s_%s', grpName, signalID));
            bpRow.Group = string(grpName);
            bpRow.LightDuration_h = meta.lightDurationHours;
            bpRow.N_Mice = 1;
            bandAll{end+1} = bpRow; %#ok<AGROW>
        end
    end

    periodTable = wavelet_vertcat_tables(periodAll);
    bandTable = wavelet_vertcat_tables(bandAll);
    writetable(periodTable, fullfile(peakFolder, sprintf('PeriodPeaks_%s.xlsx', fileStem)));

    payload = struct();
    payload.meta = meta;
    payload.groups = groups;
    payload.hsub = hsubVal;
    payload.periodTable = periodTable;
    payload.bandPower = bandTable;
    payload.figurePaths = figurePaths;
    payload.paths = struct('raw', rawFile, 'hsub', hsubFolder, 'wavelet', runFolder);
    payload.analysisNote = 'Primary wavelet on RAW; HSub is validation layer only.';

    summaryPath = write_core_summary(handoffDir, fileStem, payload);
    update_handoff_index(handoffDir, table(string(fileStem), meta.lightDurationHours, ...
        string(meta.cohort), string(summaryPath), datetime('now'), ...
        'VariableNames', {'FileStem', 'LightDuration_h', 'Cohort', 'SummaryPath', 'Created'}));

    summary = struct('runFolder', runFolder, 'handoff', summaryPath, ...
        'periodTable', periodTable, 'groups', groups);
    fprintf('Wavelet complete: %s\n', runFolder);
end

function s = sanitise_filename(strIn)
    s = char(regexprep(string(strIn), '[^\w\-]', '_'));
end
