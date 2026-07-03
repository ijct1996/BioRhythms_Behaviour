function across_render_stitched_hsub_two_panel(wtRem, wtRes, periods_hours, stitch, ...
    titleStr, outPath, theme, cfg)
%ACROSS_RENDER_STITCHED_HSUB_TWO_PANEL Stitched Removed|Residual HSub scalogram.

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1680 640]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    yLo = min(periods_hours);
    yHi = max(periods_hours);

    render_panel(nexttile(tl, 1), wtRem, periods_hours, stitch, yLo, yHi, 'Removed', theme, cfg);
    render_panel(nexttile(tl, 2), wtRes, periods_hours, stitch, yLo, yHi, 'Residual', theme, cfg);

    sgtitle(fig, titleStr, 'Interpreter', 'none', 'FontName', 'Times New Roman');
    export_figure(fig, outPath, theme, 'isScalogram', true);
    close(fig);
end

function render_panel(ax, wt, periods_hours, stitch, yLo, yHi, panelLabel, theme, cfg)
    axes(ax); %#ok<LAXES>
    pcolor(stitch.time_day, periods_hours, abs(wt));
    shading interp;
    colormap(ax, theme.scalogramColormap);
    colorbar;
    caxis auto;
    apply_plot_theme(ax, theme, 'isScalogram', true);
    set(ax, 'YTick', cfg.wavelet.scalogramYTicksHours);

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
            'Color', 'w', 'HorizontalAlignment', 'center', ...
            'FontName', 'Times New Roman', 'FontSize', 9);
    end
    hold(ax, 'off');

    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Period (hr)');
    title(ax, panelLabel);
end
