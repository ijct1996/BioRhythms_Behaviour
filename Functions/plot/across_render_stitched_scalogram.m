function across_render_stitched_scalogram(wt, periods_hours, stitch, titleStr, outPath, theme, cfg)
%ACROSS_RENDER_STITCHED_SCALOGRAM Single-panel jet scalogram with photoperiod boundaries.

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1280 640]);
    pcolor(stitch.time_day, periods_hours, abs(wt));
    shading interp;
    colormap(theme.scalogramColormap);
    colorbar;
    caxis auto;
    ax = gca;
    apply_plot_theme(ax, theme, 'isScalogram', true);
    set(ax, 'YTick', cfg.wavelet.scalogramYTicksHours);

    yLo = min(periods_hours);
    yHi = max(periods_hours);
    hold on;
    for k = 1:numel(stitch.photoBounds)
        plot([stitch.photoBounds(k) stitch.photoBounds(k)], [yLo yHi], 'w:', 'LineWidth', 2);
    end
    for k = 1:numel(stitch.lightBounds)
        plot([stitch.lightBounds(k) stitch.lightBounds(k)], [yLo yHi], 'w:', 'LineWidth', 1.5);
    end
    for k = 1:numel(stitch.segmentMid)
        text(stitch.segmentMid(k), yHi * 0.96, sprintf('L%g', stitch.photoHours(k)), ...
            'Color', 'w', 'HorizontalAlignment', 'center', ...
            'FontName', 'Times New Roman', 'FontSize', 9);
    end
    hold off;

    xlabel('Time (days)');
    ylabel('Period (hr)');
    title(titleStr, 'Interpreter', 'none');
    export_figure(fig, outPath, theme, 'isScalogram', true);
    close(fig);
end
