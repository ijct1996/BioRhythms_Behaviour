function across_plot_period_comparison(periodTable, outPath, theme, varargin)
%ACROSS_PLOT_PERIOD_COMPARISON Top peak periods across photoperiod.

    if nargin < 3, theme = plot_config('development'); end

    p = inputParser;
    addParameter(p, 'titleStr', 'Period comparison across photoperiod (exploratory)', @ischar);
    addParameter(p, 'ylabelStr', 'Dominant period (h) — rank 1', @ischar);
    addParameter(p, 'ylimRange', [0 26], @isnumeric);
    parse(p, varargin{:});

    if isempty(periodTable) || height(periodTable) == 0
        warning('across_plot_period_comparison:Empty', 'No data to plot.');
        return;
    end

    top1 = periodTable(periodTable.PeakRank == 1, :);
    if isempty(top1), return; end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 500]);
    ax = axes(fig);
    groups = unique(top1.Group);
    pal = theme.palette;
    ptSize = 36;

    hold(ax, 'on');
    for g = 1:numel(groups)
        sub = top1(top1.Group == groups(g), :);
        st = pick_collaborator_group_style(groups(g), pal);
        scatter(ax, sub.LightDuration_h, sub.PeakPeriod_hr, ptSize, ...
            'Marker', st.marker, ...
            'MarkerEdgeColor', st.markerEdgeColor, ...
            'MarkerFaceColor', st.markerFaceColor, ...
            'LineWidth', 1, ...
            'DisplayName', char(groups(g)));
    end
    hold(ax, 'off');

    xlabel(ax, 'Light duration (h)');
    ylabel(ax, p.Results.ylabelStr);
    title(ax, p.Results.titleStr);
    ylim(ax, p.Results.ylimRange);
    apply_across_xlim_padding(ax, top1.LightDuration_h);
    legend(ax, 'Location', 'best');
    apply_across_figure_theme(fig, theme);
    export_figure(fig, outPath, theme);
    close(fig);
end
