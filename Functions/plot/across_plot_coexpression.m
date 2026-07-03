function across_plot_coexpression(ratioTable, outPath, theme, varargin)
%ACROSS_PLOT_COEXPRESSION CR/UR co-expression ratios across photoperiod.

    if nargin < 3, theme = plot_config('development'); end

    p = inputParser;
    addParameter(p, 'titleStr', 'Co-expression: circadian vs ultradian band power', @ischar);
    addParameter(p, 'ylabelStr', 'CR_{20-28} / UR_{total}', @ischar);
    parse(p, varargin{:});

    if isempty(ratioTable) || height(ratioTable) == 0
        warning('across_plot_coexpression:Empty', 'No data to plot.');
        return;
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 500]);
    groups = unique(ratioTable.Group);
    pal = theme.palette;

    hold on;
    for g = 1:numel(groups)
        mask = ratioTable.Group == groups(g);
        sub = ratioTable(mask, :);
        x = unique(sub.LightDuration_h);
        mu = zeros(size(x));
        for k = 1:numel(x)
            mu(k) = mean(sub.CR_to_UR_total(sub.LightDuration_h == x(k)), 'omitnan');
        end
        st = pick_collaborator_group_style(groups(g), pal);
        plot(x, mu, '-', 'Color', st.lineColor, 'LineWidth', st.lineWidth, ...
            'Marker', st.marker, 'MarkerSize', st.markerSize, ...
            'MarkerFaceColor', st.markerFaceColor, ...
            'MarkerEdgeColor', st.markerEdgeColor, ...
            'DisplayName', char(groups(g)));
    end
    hold off;

    xlabel('Light duration (h)');
    ylabel(p.Results.ylabelStr, 'Interpreter', 'tex');
    title(p.Results.titleStr);
    legend('Location', 'best');
    apply_across_figure_theme(fig, theme);
    export_figure(fig, outPath, theme);
    close(fig);
end
