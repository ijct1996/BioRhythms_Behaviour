function outPath = wavelet_plot_scalogram(wt, periods_hours, time_day, condChangeIdx, titleStr, outPath, theme)
%WAVELET_PLOT_SCALOGRAM Jet scalogram export (development PNG / publication JPEG).

    if nargin < 7 || isempty(theme)
        theme = plot_config('development');
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 640]);
    pcolor(time_day, periods_hours, abs(wt));
    shading interp;
    colormap(theme.scalogramColormap);
    colorbar;
    caxis auto;
    ax = gca;
    apply_plot_theme(ax, theme, 'isScalogram', true);
    set(ax, 'YTick', core_defaults().wavelet.scalogramYTicksHours);

    hold on;
    for k = 1:numel(condChangeIdx)
        row = condChangeIdx(k) + 1;
        if row >= 1 && row <= numel(time_day)
            xline = time_day(row);
            plot([xline xline], [min(periods_hours) max(periods_hours)], 'w:', 'LineWidth', 1.5);
        end
    end
    hold off;

    xlabel('Time (days)');
    ylabel('Period (hr)');
    title(titleStr, 'Interpreter', 'none');

    export_figure(fig, outPath, theme, 'isScalogram', true);
    close(fig);
end
