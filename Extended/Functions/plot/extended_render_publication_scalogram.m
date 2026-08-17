function extended_render_publication_scalogram(wt, periods_hours, stitch, outPath, theme, cfg, clim)
%EXTENDED_RENDER_PUBLICATION_SCALOGRAM Single-panel stitched scalogram (Script 14).
%
%   Publication layout: Arial, jet, photoperiod labels above heatmap (not on field),
%   fixed clim, no in-image title. Figure canvas = theme.figureSizePx (default 1280×640).

    if nargin < 7 || isempty(clim) || numel(clim) ~= 2
        clim = [0, max(abs(wt(:)), [], 'omitnan')];
        if ~isfinite(clim(2)) || clim(2) <= 0
            clim(2) = 1;
        end
    end

    figW = double(theme.figureSizePx(1));
    figH = double(theme.figureSizePx(2));
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 figW figH]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on');

    pcolor(ax, stitch.time_day, periods_hours, abs(wt));
    shading(ax, 'interp');
    colormap(ax, theme.scalogramColormap);
    cb = colorbar(ax);
    cb.FontName = theme.fontName;
    cb.FontSize = 10;
    caxis(ax, clim);

    yLo = min(periods_hours);
    yHi = max(periods_hours);
    set(ax, 'YTick', cfg.wavelet.scalogramYTicksHours, ...
        'FontName', theme.fontName, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
    if ~isempty(stitch.time_day)
        set(ax, 'XLim', [min(stitch.time_day), max(stitch.time_day)]);
    end
    set(ax, 'YLim', [yLo, yHi]);

    for k = 1:numel(stitch.photoBounds)
        x = stitch.photoBounds(k);
        plot(ax, [x x], [yLo yHi], 'w:', 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    for k = 1:numel(stitch.lightBounds)
        x = stitch.lightBounds(k);
        plot(ax, [x x], [yLo yHi], 'w:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    for k = 1:numel(stitch.segmentMid)
        text(ax, stitch.segmentMid(k), yHi, extended_script14_pp_label_(stitch.photoHours(k)), ...
            'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontName', theme.fontName, 'FontSize', 9, 'Clipping', 'off', ...
            'Interpreter', 'none', 'HandleVisibility', 'off');
    end
    hold(ax, 'off');

    xlabel(ax, 'Time (days)', 'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 12);
    ylabel(ax, 'Period (hr)', 'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 12);

    export_figure(fig, outPath, theme, 'isScalogram', true);
    close(fig);
end

function lbl = extended_script14_pp_label_(photoH)
    photoH = round(double(photoH));
    if photoH == 24
        lbl = 'LL';
    else
        lbl = sprintf('L%d', photoH);
    end
end
