function stitch = extended_script14_stitch_raw_group(entries, grpName, cfg)
%EXTENDED_SCRIPT14_STITCH_RAW_GROUP Stitched group-mean RAW across photoperiods.
%   Adapted from Core across_plot_stitched_scalograms (Script 3); Extended-only.

    stitch = script14_init_stitch_();
    dayOffset = 0;

    for i = 1:numel(entries)
        grp = script14_find_group_(entries(i).groups, grpName);
        if isempty(grp)
            warning('extended_script14:MissingGroup', ...
                'Group "%s" not defined in %s; segment skipped.', grpName, entries(i).fileStem);
            continue;
        end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'raw')
            error('extended_script14:NoRawPath', ...
                'Handoff summary for %s lacks paths.raw (re-run Script 2).', entries(i).fileStem);
        end
        rawPath = entries(i).paths.raw;
        if ~isfile(rawPath)
            error('extended_script14:MissingRaw', 'RAW file not found: %s', rawPath);
        end

        [tbl, meta] = read_behav_excel(rawPath);
        grp = resolve_groups_by_names(meta.varNames, grp);
        if isempty(grp.colIdx)
            continue;
        end

        segSignal = mean(wavelet_stack_signals(tbl, grp.colIdx), 2);
        n = numel(segSignal);
        segDay = ((0:n - 1)' * cfg.samplingMinutes) / (60 * 24) + dayOffset;
        lightVec = tbl{:, meta.lightIdx};
        condChangeIdx = find(diff(lightVec) ~= 0);
        lightBounds = [];
        for k = 1:numel(condChangeIdx)
            row = condChangeIdx(k) + 1;
            if row >= 1 && row <= n
                lightBounds(end + 1) = segDay(row); %#ok<AGROW>
            end
        end

        stitch = script14_append_segment_(stitch, segSignal, segDay, dayOffset, ...
            entries(i).lightHours, lightBounds, cfg);
        dayOffset = stitch.nextDayOffset;
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
