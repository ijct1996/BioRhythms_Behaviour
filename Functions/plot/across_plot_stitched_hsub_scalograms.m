function figurePaths = across_plot_stitched_hsub_scalograms(entries, figFolder, cohortTag, theme)
%ACROSS_PLOT_STITCHED_HSUB_SCALOGRAMS Stitched SEL_P360 Removed|Residual HSub scalograms.
%
%   Reads Script 1 TimeSeries exports via paths.hsub in each handoff summary.

    if nargin < 4 || isempty(theme)
        theme = plot_config('development');
    end
    theme = plot_theme_ensure_scalogram(theme);
    cfg = core_defaults();
    armKey = cfg.hsub.scalogramArm;
    armLabel = cfg.hsub.scalogramLabel;

    groupNames = across_collect_group_names(entries);
    if isempty(groupNames)
        warning('across_plot_stitched_hsub_scalograms:NoGroups', 'No condition groups in handoff.');
        figurePaths = {};
        return;
    end

    stitchedFolder = fullfile(figFolder, 'Stitched_Scalograms_HSub');
    ensure_dir(stitchedFolder);
    ext = theme.scalogram.format;
    figurePaths = {};

    for g = 1:numel(groupNames)
        grpName = groupNames{g};
        [stitchRem, stitchRes] = across_stitch_hsub_group(entries, grpName, armKey, cfg);
        if isempty(stitchRem.signal) || isempty(stitchRes.signal)
            warning('across_plot_stitched_hsub_scalograms:EmptySignal', ...
                'Skipping group "%s": no valid HSub segments.', grpName);
            continue;
        end

        FB = wavelet_make_filterbank(numel(stitchRem.signal), cfg.samplingMinutes, cfg);
        [wtRem, periods_hours, ~] = wavelet_compute_cwt(stitchRem.signal, FB);
        [wtRes, ~, ~] = wavelet_compute_cwt(stitchRes.signal, FB);

        safeName = sanitise_filename(grpName);
        outFile = fullfile(stitchedFolder, sprintf('Stitched_HSub_Removed_Residual_%s_%s_%s.%s', ...
            armLabel, safeName, cohortTag, ext));
        titleStr = sprintf('HSub %s | Stitched average | %s | %s | L%g–L%g h', ...
            armLabel, grpName, cohortTag, min(stitchRem.photoHours), max(stitchRem.photoHours));

        across_render_stitched_hsub_two_panel(wtRem, wtRes, periods_hours, stitchRem, ...
            titleStr, outFile, theme, cfg);
        figurePaths{end + 1} = outFile; %#ok<AGROW>
        fprintf('  Stitched HSub scalogram: "%s" (%d photoperiods) → %s\n', ...
            grpName, numel(stitchRem.photoHours), outFile);
    end
end

function [stitchRem, stitchRes] = across_stitch_hsub_group(entries, grpName, armKey, cfg)
    stitchRem = init_stitch();
    stitchRes = init_stitch();
    dayOffset = 0;

    for i = 1:numel(entries)
        grp = find_group(entries(i).groups, grpName);
        if isempty(grp)
            warning('across_stitch_hsub:MissingGroup', ...
                'Group "%s" not defined in %s; segment skipped.', grpName, entries(i).fileStem);
            continue;
        end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'hsub')
            error('across_stitch_hsub:NoHSubPath', ...
                'Handoff summary for %s lacks paths.hsub (re-run Script 2).', entries(i).fileStem);
        end

        hsubRun = fullfile(entries(i).paths.hsub, entries(i).fileStem);
        [avgRem, avgRes, segDay, lightBounds] = hsub_load_group_arm_average( ...
            hsubRun, entries(i).fileStem, grp, armKey);
        if isempty(avgRem)
            continue;
        end

        segDay = segDay + dayOffset;
        stitchRem = append_segment(stitchRem, avgRem, segDay, dayOffset, entries(i).lightHours, lightBounds);
        stitchRes = append_segment(stitchRes, avgRes, segDay, dayOffset, entries(i).lightHours, lightBounds);
        dayOffset = stitchRem.nextDayOffset;
    end
end

function [avgRem, avgRes, time_day, lightBounds] = hsub_load_group_arm_average( ...
    hsubRunFolder, fileStem, grp, armKey)

    avgRem = [];
    avgRes = [];
    time_day = [];
    lightBounds = [];

    remFile = fullfile(hsubRunFolder, 'TimeSeries', 'Removed', ...
        sprintf('Removed_Selective_%s_%s.xlsx', armKey, fileStem));
    resFile = fullfile(hsubRunFolder, 'TimeSeries', 'Residual', ...
        sprintf('Residual_Selective_%s_%s.xlsx', armKey, fileStem));
    if ~isfile(remFile) || ~isfile(resFile)
        warning('hsub_load_group_arm_average:MissingFile', ...
            'HSub timeseries missing for %s (re-run Script 1).', fileStem);
        return;
    end

    remTbl = readtable(remFile);
    resTbl = readtable(resFile);
    colNames = grp.colNames;
    if isstring(colNames), colNames = cellstr(colNames); end

    remCols = [];
    resCols = [];
    for c = 1:numel(colNames)
        remVar = matlab.lang.makeValidName(['Removed_Selective_' colNames{c}]);
        resVar = matlab.lang.makeValidName(['Residual_Selective_' colNames{c}]);
        if ismember(remVar, remTbl.Properties.VariableNames)
            remCols(:, end + 1) = remTbl.(remVar); %#ok<AGROW>
        end
        if ismember(resVar, resTbl.Properties.VariableNames)
            resCols(:, end + 1) = resTbl.(resVar); %#ok<AGROW>
        end
    end
    if isempty(remCols) || isempty(resCols)
        return;
    end

    avgRem = mean(remCols, 2, 'omitnan');
    avgRes = mean(resCols, 2, 'omitnan');
    n = numel(avgRem);
    cfg = core_defaults();
    time_day = ((0:n - 1)' * cfg.samplingMinutes) / (60 * 24);

    lightVar = remTbl.Properties.VariableNames{end};
    if contains(lower(lightVar), 'light')
        lightVec = remTbl.(lightVar);
        condChangeIdx = find(diff(lightVec) ~= 0);
        lightBounds = [];
        for k = 1:numel(condChangeIdx)
            row = condChangeIdx(k) + 1;
            if row >= 1 && row <= n
                lightBounds(end + 1) = time_day(row); %#ok<AGROW>
            end
        end
    else
        lightBounds = [];
    end
end

function names = across_collect_group_names(entries)
    names = {};
    for i = 1:numel(entries)
        if ~isfield(entries(i), 'groups') || isempty(entries(i).groups)
            continue;
        end
        for g = 1:numel(entries(i).groups)
            names{end + 1} = entries(i).groups(g).name; %#ok<AGROW>
        end
    end
    names = unique(names, 'stable');
end

function stitch = init_stitch()
    stitch = struct('signal', [], 'time_day', [], 'photoBounds', [], ...
        'lightBounds', [], 'photoHours', [], 'segmentMid', [], 'nextDayOffset', 0);
end

function stitch = append_segment(stitch, segSignal, segDay, dayOffset, lightHours, lightBounds)
    cfg = core_defaults();
    TsDay = cfg.samplingMinutes / (60 * 24);
    segStart = dayOffset;

    if ~isempty(stitch.signal)
        stitch.photoBounds(end + 1) = dayOffset; %#ok<AGROW>
    end
    if ~isempty(lightBounds)
        stitch.lightBounds = [stitch.lightBounds, lightBounds(:)']; %#ok<AGROW>
    end

    stitch.signal = [stitch.signal; segSignal(:)]; %#ok<AGROW>
    stitch.time_day = [stitch.time_day; segDay(:)]; %#ok<AGROW>
    stitch.photoHours(end + 1) = lightHours; %#ok<AGROW>
    stitch.segmentMid(end + 1) = (segStart + segDay(end)) / 2; %#ok<AGROW>
    stitch.nextDayOffset = segDay(end) + TsDay;
end

function grp = find_group(groups, grpName)
    grp = [];
    for k = 1:numel(groups)
        if strcmp(groups(k).name, grpName)
            grp = groups(k);
            return;
        end
    end
end

function s = sanitise_filename(strIn)
    s = char(regexprep(string(strIn), '[^\w\-]', '_'));
end
