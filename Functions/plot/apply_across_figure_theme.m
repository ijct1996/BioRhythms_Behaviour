function apply_across_figure_theme(fig, theme)
%APPLY_ACROSS_FIGURE_THEME Collaborator typography for Script 3 line/scatter figures.

    if nargin < 2 || isempty(theme)
        theme = plot_config('development');
    end

    axList = findall(fig, 'Type', 'axes');
    for k = 1:numel(axList)
        ax = axList(k);
        apply_plot_theme(ax, theme);
        if strcmp(theme.mode, 'publication')
            style_text(get(ax, 'XLabel'), theme);
            style_text(get(ax, 'YLabel'), theme);
            style_text(get(ax, 'Title'), theme);
        end
    end

    if strcmp(theme.mode, 'publication')
        leg = findobj(fig, 'Type', 'Legend');
        if ~isempty(leg)
            set(leg, 'FontName', theme.axesFontName, ...
                'FontSize', theme.axesFontSize, ...
                'FontWeight', theme.axesFontWeight);
        end
    end
end

function style_text(h, theme)
    if ~isempty(h) && isgraphics(h)
        set(h, 'FontName', theme.axesFontName, ...
            'FontSize', theme.axesFontSize, ...
            'FontWeight', theme.axesFontWeight);
    end
end
