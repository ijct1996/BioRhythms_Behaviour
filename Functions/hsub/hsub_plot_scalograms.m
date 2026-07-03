function figurePaths = hsub_plot_scalograms(tbl, meta, groups, residualMat, removedMat, ...
    mouseIdx, anchorOK, plotMode, saveIndividual, fileStem, runFolder)
%HSUB_PLOT_SCALOGRAMS Two-panel Removed|Residual group-average HSub scalograms.
%
%   Uses SEL_P360 (Selective Min360) signals. Anchor-rejected mice are excluded
%   from group averages. Optional individual two-panel figures at 150 DPI PNG.

    cfg = core_defaults();
    theme = plot_config(plotMode);
    theme = plot_theme_ensure_scalogram(theme);

    figurePaths = {};
    if ~cfg.hsub.saveScalograms
        return;
    end

    time_hours = tbl{:, meta.timeIdx};
    time_min = time_hours * 60;
    time_day = time_min / (60 * 24);
    TsMinutes = median(diff(time_min), 'omitnan');
    lightVec = tbl{:, meta.lightIdx};
    condChangeIdx = find(diff(lightVec) ~= 0);

    scaloRoot = fullfile(runFolder, 'Figures_Scalograms');
    ensure_dir(scaloRoot);
    ext = theme.scalogram.format;
    armLabel = cfg.hsub.scalogramLabel;
    photoTag = photoperiod_tag(meta);

    FB = wavelet_make_filterbank(height(tbl), TsMinutes, cfg);

    for g = 1:numel(groups)
        grp = resolve_groups_by_names(meta.varNames, groups(g));
        if isempty(grp.colIdx)
            warning('hsub_plot_scalograms:EmptyGroup', ...
                'Group "%s" has no columns in %s; skipped.', groups(g).name, fileStem);
            continue;
        end

        grpFolder = fullfile(scaloRoot, ['Group_' sanitise_filename(grp.name)]);
        ensure_dir(grpFolder);

        colMask = false(1, numel(mouseIdx));
        for c = 1:numel(mouseIdx)
            if ~anchorOK(c)
                continue;
            end
            if any(mouseIdx(c) == grp.colIdx)
                colMask(c) = true;
            end
        end
        nOk = sum(colMask);
        if nOk == 0
            warning('hsub_plot_scalograms:NoAnchorOK', ...
                'Group "%s" has no AnchorOK mice in %s; skipped.', grp.name, fileStem);
            continue;
        end

        avgRemoved = mean(removedMat(:, colMask), 2);
        avgResidual = mean(residualMat(:, colMask), 2);
        [wtRem, periods_hours, ~] = wavelet_compute_cwt(avgRemoved, FB);
        [wtRes, ~, ~] = wavelet_compute_cwt(avgResidual, FB);

        safeName = sanitise_filename(grp.name);
        outFile = fullfile(grpFolder, sprintf('HSub_Removed_Residual_%s_%s_%s.%s', ...
            armLabel, safeName, photoTag, ext));
        titleStr = sprintf('HSub %s | Average | %s | n=%d AnchorOK | L=%g h', ...
            armLabel, grp.name, nOk, meta.lightDurationHours);

        hsub_plot_two_panel_scalogram(wtRem, wtRes, periods_hours, time_day, ...
            condChangeIdx, titleStr, outFile, theme, cfg);
        figurePaths{end + 1} = outFile; %#ok<AGROW>
        fprintf('  HSub scalogram: "%s" (n=%d) → %s\n', grp.name, nOk, outFile);

        if saveIndividual
            indivFolder = fullfile(grpFolder, 'Individuals');
            ensure_dir(indivFolder);
            qcCfg = cfg.plot.hsubIndividual;
            for c = find(colMask)
                mouseName = meta.varNames{mouseIdx(c)};
                safeMouse = sanitise_filename(mouseName);
                sigRem = removedMat(:, c);
                sigRes = residualMat(:, c);
                [wtRemI, ~, ~] = wavelet_compute_cwt(sigRem, FB);
                [wtResI, ~, ~] = wavelet_compute_cwt(sigRes, FB);
                outIndiv = fullfile(indivFolder, sprintf('HSub_Removed_Residual_%s_%s_%s.%s', ...
                    armLabel, safeName, safeMouse, qcCfg.format));
                titleIndiv = sprintf('HSub %s | %s | %s | L=%g h', ...
                    armLabel, grp.name, mouseName, meta.lightDurationHours);
                hsub_plot_two_panel_scalogram(wtRemI, wtResI, periods_hours, time_day, ...
                    condChangeIdx, titleIndiv, outIndiv, theme, cfg, ...
                    'dpiOverride', qcCfg.dpi, 'formatOverride', qcCfg.format);
            end
        end
    end
end

function hsub_plot_two_panel_scalogram(wtRem, wtRes, periods_hours, time_day, ...
    condChangeIdx, titleStr, outPath, theme, cfg, varargin)

    p = inputParser;
    addParameter(p, 'dpiOverride', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'formatOverride', '', @ischar);
    parse(p, varargin{:});

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1600 640]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    yLo = min(periods_hours);
    yHi = max(periods_hours);

    ax1 = nexttile(tl, 1);
    plot_scalogram_panel(ax1, wtRem, periods_hours, time_day, condChangeIdx, yLo, yHi, ...
        'Removed', theme, cfg);
    ax2 = nexttile(tl, 2);
    plot_scalogram_panel(ax2, wtRes, periods_hours, time_day, condChangeIdx, yLo, yHi, ...
        'Residual', theme, cfg);

    sgtitle(fig, titleStr, 'Interpreter', 'none', 'FontName', 'Times New Roman');
    export_figure(fig, outPath, theme, 'isScalogram', true, ...
        'dpiOverride', p.Results.dpiOverride, 'formatOverride', p.Results.formatOverride);
    close(fig);
end

function plot_scalogram_panel(ax, wt, periods_hours, time_day, condChangeIdx, yLo, yHi, ...
    panelLabel, theme, cfg)

    axes(ax); %#ok<LAXES>
    pcolor(time_day, periods_hours, abs(wt));
    shading interp;
    colormap(ax, theme.scalogramColormap);
    colorbar;
    caxis auto;
    apply_plot_theme(ax, theme, 'isScalogram', true);
    set(ax, 'YTick', cfg.wavelet.scalogramYTicksHours);

    hold(ax, 'on');
    for k = 1:numel(condChangeIdx)
        row = condChangeIdx(k) + 1;
        if row >= 1 && row <= numel(time_day)
            xline = time_day(row);
            plot(ax, [xline xline], [yLo yHi], 'w:', 'LineWidth', 1.5);
        end
    end
    hold(ax, 'off');

    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Period (hr)');
    title(ax, panelLabel);
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
