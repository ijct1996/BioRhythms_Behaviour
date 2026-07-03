function summary = hsub_run_file(inputFile, outputFolder, varargin)
%HSUB_RUN_FILE Run harmonic subtraction on one Pipeline_Input xlsx.
%
%   summary = hsub_run_file(inputFile, outputFolder)
%   summary = hsub_run_file(inputFile, outputFolder, 'Interactive', false, ...
%       'MouseIdx', mouseIdx, 'Meta', meta)

    p = inputParser;
    addParameter(p, 'Interactive', true, @islogical);
    addParameter(p, 'MouseIdx', [], @isnumeric);
    addParameter(p, 'MouseNames', {}, @(x) iscell(x) || isstring(x));
    addParameter(p, 'Meta', struct(), @isstruct);
    addParameter(p, 'PlotMode', 'development', @ischar);
    addParameter(p, 'Groups', [], @(x) isempty(x) || isstruct(x));
    addParameter(p, 'SaveIndividualScalograms', false, @islogical);
    parse(p, varargin{:});

    [tbl, meta] = read_behav_excel(inputFile);
    if ~isempty(fieldnames(p.Results.Meta))
        meta = p.Results.Meta;
    end

    opts = hsub_get_opts();
    [timeMinutesAll, TsMinutes] = hsub_infer_time(tbl{:, meta.timeIdx}, meta.timeName);
    opts.TsMinutes = TsMinutes;

    durationHours = (max(timeMinutesAll) - min(timeMinutesAll)) / 60;
    baselineWinHours = min(max(72, 0.5 * durationHours), 168);
    opts.baselineWinSamples = max(20, round((baselineWinHours * 60) / TsMinutes));
    opts.blockLenSamples = max(4, round((opts.blockLenHours * 60) / TsMinutes));

    mouseIdx = p.Results.MouseIdx;
    if isempty(mouseIdx)
        mouseIdx = meta.mouseIdx;
    end

    mouseNamesArg = p.Results.MouseNames;
    if isstring(mouseNamesArg), mouseNamesArg = cellstr(mouseNamesArg); end

    if ~isempty(mouseNamesArg)
        mouseIdx = resolve_column_indices(meta.varNames, mouseNamesArg);
        fprintf('Reusing %d mouse columns from batch template.\n', numel(mouseIdx));
    elseif p.Results.Interactive
        defaultData = mouseIdx;
        [dataIdx, ok] = listdlg('PromptString', 'Select mouse columns:', ...
            'SelectionMode', 'multiple', 'ListString', meta.varNames, ...
            'InitialValue', defaultData);
        if ~ok || isempty(dataIdx), error('hsub_run_file:NoColumns', 'No columns selected.'); end
        mouseIdx = dataIdx;
    end

    mouseNames = meta.varNames(mouseIdx);

    [~, fileName, ~] = fileparts(inputFile);
    runFolder = fullfile(outputFolder, fileName);
    reportsFolder = fullfile(runFolder, 'Reports');
    tsResFolder = fullfile(runFolder, 'TimeSeries', 'Residual');
    tsRemFolder = fullfile(runFolder, 'TimeSeries', 'Removed');
    figAnchorFolder = fullfile(runFolder, 'Figures_Anchor');

    ensure_dir(reportsFolder);
    ensure_dir(tsResFolder);
    ensure_dir(tsRemFolder);
    ensure_dir(figAnchorFolder);

    lightVec = tbl{:, meta.lightIdx};
    timeOut = tbl{:, meta.timeIdx};
    lightOut = lightVec;
    LightDur_h = meta.lightDurationHours;

    nMice = numel(mouseIdx);
    anchorRows = cell(nMice, 8);
    recRows = cell(nMice, 7);
    results = cell(nMice, 1);

    keys = arrayfun(@(mp) sprintf('Min%d', mp), opts.minPeriodsToRun, 'UniformOutput', false);
    resArrays = struct();
    for ki = 1:numel(keys)
        resArrays.(keys{ki}).sel = NaN(height(tbl), nMice);
        resArrays.(keys{ki}).full = NaN(height(tbl), nMice);
        resArrays.(keys{ki}).remSel = NaN(height(tbl), nMice);
        resArrays.(keys{ki}).remFull = NaN(height(tbl), nMice);
    end

    for c = 1:nMice
        colIdx = mouseIdx(c);
        colName = meta.varNames{colIdx};
        fprintf('HSub %d/%d: %s\n', c, nMice, colName);

        x0 = tbl{:, colIdx};
        if ~isnumeric(x0), x0 = str2double(string(x0)); end

        mouseOut = hsub_process_mouse(x0, timeMinutesAll, opts, LightDur_h);
        mouseOut.colName = colName;
        rec = hsub_build_recommendation(mouseOut, opts);
        results{c} = struct('mouse', mouseOut, 'rec', rec);

        anchorRows(c, :) = {colName, mouseOut.anchorOK, mouseOut.P0_h, ...
            mouseOut.maxDeltaR2, mouseOut.pBlock, mouseOut.cyclesAtP0, mouseOut.note, rec.code};

        recRows(c, :) = {colName, rec.code, rec.mode, rec.minKey, rec.anchorOK, rec.reason, mouseOut.P0_h};

        for ki = 1:numel(keys)
            k = keys{ki};
            if isfield(mouseOut.residuals, k)
                resArrays.(k).sel(:, c) = mouseOut.residuals.(k).selective;
                resArrays.(k).full(:, c) = mouseOut.residuals.(k).full;
                resArrays.(k).remSel(:, c) = mouseOut.residuals.(k).removedSel;
                resArrays.(k).remFull(:, c) = mouseOut.residuals.(k).removedFull;
            end
        end
    end

    % Write four residual Excel arms + recommended
    for ki = 1:numel(keys)
        k = keys{ki};
        hsub_write_timeseries(tsResFolder, sprintf('Residual_Selective_%s_%s.xlsx', k, fileName), ...
            meta, timeOut, lightOut, mouseIdx, resArrays.(k).sel, 'Residual_Selective');
        hsub_write_timeseries(tsResFolder, sprintf('Residual_FullLadder_%s_%s.xlsx', k, fileName), ...
            meta, timeOut, lightOut, mouseIdx, resArrays.(k).full, 'Residual_FullLadder');
        hsub_write_timeseries(tsRemFolder, sprintf('Removed_Selective_%s_%s.xlsx', k, fileName), ...
            meta, timeOut, lightOut, mouseIdx, resArrays.(k).remSel, 'Removed_Selective');
        hsub_write_timeseries(tsRemFolder, sprintf('Removed_FullLadder_%s_%s.xlsx', k, fileName), ...
            meta, timeOut, lightOut, mouseIdx, resArrays.(k).remFull, 'Removed_FullLadder');
    end

    % Recommended residual per mouse (SEL_P360 default path)
    recMatrix = NaN(height(tbl), nMice);
    for c = 1:nMice
        r = results{c}.rec;
        k = r.minKey;
        if strcmp(r.mode, 'FullLadder')
            recMatrix(:, c) = resArrays.(k).full(:, c);
        elseif r.anchorOK
            recMatrix(:, c) = resArrays.(k).sel(:, c);
        else
            recMatrix(:, c) = tbl{:, mouseIdx(c)};
        end
    end
    hsub_write_timeseries(tsResFolder, sprintf('Residual_Recommended_%s.xlsx', fileName), ...
        meta, timeOut, lightOut, mouseIdx, recMatrix, 'Residual_Recommended');

    anchorTable = cell2table(anchorRows, 'VariableNames', ...
        {'Column', 'AnchorOK', 'Period_hours', 'MaxDeltaR2', 'pBlock', 'Cycles', 'Notes', 'Recommendation'});
    recTable = cell2table(recRows, 'VariableNames', ...
        {'Column', 'Code', 'Mode', 'MinKey', 'AnchorOK', 'Reason', 'Period_hours'});

    summaryXlsx = fullfile(reportsFolder, 'HarmonicRemoval_Summary.xlsx');
    if isfile(summaryXlsx), delete(summaryXlsx); end
    writetable(anchorTable, summaryXlsx, 'Sheet', 'Anchor_Report');
    writetable(recTable, summaryXlsx, 'Sheet', 'Recommendation');

    groups = p.Results.Groups;
    if isempty(groups)
        groups = groups_to_template(group_assignment_dialog( ...
            meta.varNames, mouseIdx, meta, 'Interactive', false), meta.varNames);
    end

    anchorOK = false(nMice, 1);
    for c = 1:nMice
        anchorOK(c) = results{c}.mouse.anchorOK;
    end

    hsubArm = core_defaults().hsub.scalogramArm;
    if isfield(resArrays, hsubArm)
        hsubFigPaths = hsub_plot_scalograms(tbl, meta, groups, ...
            resArrays.(hsubArm).sel, resArrays.(hsubArm).remSel, ...
            mouseIdx, anchorOK, p.Results.PlotMode, p.Results.SaveIndividualScalograms, ...
            fileName, runFolder);
    else
        hsubFigPaths = {};
    end

    summary = struct();
    summary.inputFile = inputFile;
    summary.outputFolder = runFolder;
    summary.meta = meta;
    summary.mouseIdx = mouseIdx;
    summary.mouseNames = mouseNames;
    summary.groups = groups;
    summary.results = results;
    summary.summaryXlsx = summaryXlsx;
    summary.hsubFigurePaths = hsubFigPaths;
    fprintf('HSub complete: %s\n', runFolder);
end

function hsub_write_timeseries(outFolder, outName, meta, timeOut, lightOut, mouseIdx, dataMatrix, prefix)
    headers = [{'Time'}, strcat(prefix, '_', meta.varNames(mouseIdx)), {'Light_duration_h'}];
    T = array2table([timeOut, dataMatrix, lightOut], 'VariableNames', matlab.lang.makeValidName(headers));
    writetable(T, fullfile(outFolder, outName));
end
