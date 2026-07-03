function [figurePaths, periodRows, bandRows] = wavelet_plot_group_averages( ...
    tbl, meta, groups, FB, time_day, condChangeIdx, scaloRoot, theme, cfg)
%WAVELET_PLOT_GROUP_AVERAGES One mean scalogram per user-defined condition group.
%
%   Saved to Scalograms/Group_{name}/Scalogram_Average_{name}_L{photoperiod}.{ext}

    figurePaths = {};
    periodRows = {};
    bandRows = {};
    ext = theme.scalogram.format;
    photoTag = photoperiod_tag(meta);

    for g = 1:numel(groups)
        grpName = groups(g).name;
        colIdx = groups(g).colIdx;
        if isempty(colIdx)
            continue;
        end

        grpFolder = fullfile(scaloRoot, ['Group_' sanitise_filename(grpName)]);
        ensure_dir(grpFolder);

        avgSignal = mean(wavelet_stack_signals(tbl, colIdx), 2);
        [wt, periods_hours, ~] = wavelet_compute_cwt(avgSignal, FB);
        powerSpec = mean(abs(wt).^2, 2);

        safeName = sanitise_filename(grpName);
        outFile = fullfile(grpFolder, sprintf('Scalogram_Average_%s_%s.%s', safeName, photoTag, ext));
        titleStr = sprintf('RAW | Average | %s | n=%d | L=%g h', ...
            grpName, numel(colIdx), meta.lightDurationHours);
        wavelet_plot_scalogram(wt, periods_hours, time_day, condChangeIdx, ...
            titleStr, outFile, theme);
        figurePaths{end + 1} = outFile; %#ok<AGROW>

        signalLabel = sprintf('Average_%s', safeName);
        peaks = wavelet_find_period_peaks(powerSpec, periods_hours, cfg.wavelet.topNPeaks);
        if height(peaks) > 0
            peaks.SignalID = repmat(string(signalLabel), height(peaks), 1);
            peaks.Group = repmat(string(grpName), height(peaks), 1);
            peaks.LightDuration_h = repmat(meta.lightDurationHours, height(peaks), 1);
            peaks.N_Mice = repmat(numel(colIdx), height(peaks), 1);
            periodRows{end + 1} = peaks; %#ok<AGROW>
        end

        bp = wavelet_band_power(powerSpec, periods_hours, cfg);
        bpRow = struct2table(bp);
        bpRow.SignalID = string(signalLabel);
        bpRow.Group = string(grpName);
        bpRow.LightDuration_h = meta.lightDurationHours;
        bpRow.N_Mice = numel(colIdx);
        bandRows{end + 1} = bpRow; %#ok<AGROW>

        fprintf('  Group average: "%s" (%d mice) → %s\n', grpName, numel(colIdx), outFile);
    end
end

function tag = photoperiod_tag(meta)
    if isfield(meta, 'lightDurationHours') && isfinite(meta.lightDurationHours)
        tag = sprintf('L%g', meta.lightDurationHours);
    else
        tag = meta.fileStem;
    end
end

function s = sanitise_filename(strIn)
    s = char(regexprep(string(strIn), '[^\w\-]', '_'));
end
