function figurePaths = across_plot_stitched_scalograms(entries, figFolder, cohortTag, theme)
%ACROSS_PLOT_STITCHED_SCALOGRAMS Group-average RAW scalograms stitched across photoperiod.

    if nargin < 4 || isempty(theme)
        theme = plot_config('development');
    end
    theme = plot_theme_ensure_scalogram(theme);
    cfg = core_defaults();

    groupNames = across_collect_group_names(entries);
    if isempty(groupNames)
        warning('across_plot_stitched_scalograms:NoGroups', 'No condition groups in handoff.');
        figurePaths = {};
        return;
    end

    stitchedFolder = fullfile(figFolder, 'Stitched_Scalograms');
    ensure_dir(stitchedFolder);
    ext = theme.scalogram.format;
    figurePaths = {};

    for g = 1:numel(groupNames)
        grpName = groupNames{g};
        stitch = across_stitch_raw_group(entries, grpName, cfg);
        if isempty(stitch.signal)
            warning('across_plot_stitched_scalograms:EmptySignal', ...
                'Skipping group "%s": no valid photoperiod segments.', grpName);
            continue;
        end

        FB = wavelet_make_filterbank(numel(stitch.signal), cfg.samplingMinutes, cfg);
        [wt, periods_hours, ~] = wavelet_compute_cwt(stitch.signal, FB);

        safeName = sanitise_filename(grpName);
        outFile = fullfile(stitchedFolder, sprintf('Stitched_Scalogram_Average_%s_%s.%s', ...
            safeName, cohortTag, ext));
        titleStr = sprintf('RAW | Stitched average | %s | %s | L%g–L%g h', ...
            grpName, cohortTag, min(stitch.photoHours), max(stitch.photoHours));

        across_render_stitched_scalogram(wt, periods_hours, stitch, titleStr, outFile, theme, cfg);
        figurePaths{end + 1} = outFile; %#ok<AGROW>
        fprintf('  Stitched RAW scalogram: "%s" (%d photoperiods) → %s\n', ...
            grpName, numel(stitch.photoHours), outFile);
    end
end

function stitch = across_stitch_raw_group(entries, grpName, cfg)
    stitch = init_stitch();
    dayOffset = 0;

    for i = 1:numel(entries)
        grp = find_group(entries(i).groups, grpName);
        if isempty(grp)
            warning('across_stitch:MissingGroup', ...
                'Group "%s" not defined in %s; segment skipped.', grpName, entries(i).fileStem);
            continue;
        end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'raw')
            error('across_stitch:NoRawPath', ...
                'Handoff summary for %s lacks paths.raw (re-run Script 2).', entries(i).fileStem);
        end
        rawPath = entries(i).paths.raw;
        if ~isfile(rawPath)
            error('across_stitch:MissingRaw', 'RAW file not found: %s', rawPath);
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

        stitch = append_segment(stitch, segSignal, segDay, dayOffset, entries(i).lightHours, lightBounds);
        dayOffset = stitch.nextDayOffset;
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
