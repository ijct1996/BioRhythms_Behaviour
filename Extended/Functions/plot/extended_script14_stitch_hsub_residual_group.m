function stitch = extended_script14_stitch_hsub_residual_group(entries, grpName, armKey, cfg)
%EXTENDED_SCRIPT14_STITCH_HSUB_RESIDUAL_GROUP Stitched SEL_P360 residual group mean.
%   Adapted from Core across_plot_stitched_hsub_scalograms (Script 3); residual only.

    stitch = script14_init_stitch_();
    dayOffset = 0;

    for i = 1:numel(entries)
        grp = script14_find_group_(entries(i).groups, grpName);
        if isempty(grp)
            warning('extended_script14:MissingGroup', ...
                'Group "%s" not defined in %s; segment skipped.', grpName, entries(i).fileStem);
            continue;
        end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'hsub')
            error('extended_script14:NoHSubPath', ...
                'Handoff summary for %s lacks paths.hsub (re-run Script 2).', entries(i).fileStem);
        end

        hsubRun = fullfile(entries(i).paths.hsub, entries(i).fileStem);
        [avgRes, segDay, lightBounds] = script14_hsub_load_group_residual_( ...
            hsubRun, entries(i).fileStem, grp, armKey, cfg);
        if isempty(avgRes)
            continue;
        end

        segDay = segDay + dayOffset;
        stitch = script14_append_segment_(stitch, avgRes, segDay, dayOffset, ...
            entries(i).lightHours, lightBounds, cfg);
        dayOffset = stitch.nextDayOffset;
    end
end

function [avgRes, time_day, lightBounds] = script14_hsub_load_group_residual_( ...
    hsubRunFolder, fileStem, grp, armKey, cfg)

    avgRes = [];
    time_day = [];
    lightBounds = [];

    resFile = fullfile(hsubRunFolder, 'TimeSeries', 'Residual', ...
        sprintf('Residual_Selective_%s_%s.xlsx', armKey, fileStem));
    if ~isfile(resFile)
        warning('extended_script14:MissingHSubResidual', ...
            'HSub residual timeseries missing for %s (re-run Script 1).', fileStem);
        return;
    end

    resTbl = readtable(resFile);
    colNames = grp.colNames;
    if isstring(colNames), colNames = cellstr(colNames); end

    resCols = [];
    for c = 1:numel(colNames)
        resVar = matlab.lang.makeValidName(['Residual_Selective_' colNames{c}]);
        if ismember(resVar, resTbl.Properties.VariableNames)
            resCols(:, end + 1) = resTbl.(resVar); %#ok<AGROW>
        end
    end
    if isempty(resCols)
        return;
    end

    avgRes = mean(resCols, 2, 'omitnan');
    n = numel(avgRes);
    time_day = ((0:n - 1)' * cfg.samplingMinutes) / (60 * 24);

    lightVar = resTbl.Properties.VariableNames{end};
    if contains(lower(lightVar), 'light')
        lightVec = resTbl.(lightVar);
        condChangeIdx = find(diff(lightVec) ~= 0);
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

function stitch = script14_init_stitch_()
    stitch = struct('signal', [], 'time_day', [], 'photoBounds', [], ...
        'lightBounds', [], 'photoHours', [], 'segmentMid', [], 'nextDayOffset', 0);
end

function stitch = script14_append_segment_(stitch, segSignal, segDay, dayOffset, lightHours, lightBounds, cfg)
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
    TsDay = cfg.samplingMinutes / (60 * 24);
    stitch.nextDayOffset = segDay(end) + TsDay;
end

function grp = script14_find_group_(groups, grpName)
    grp = [];
    for k = 1:numel(groups)
        if strcmp(groups(k).name, grpName)
            grp = groups(k);
            return;
        end
    end
end
