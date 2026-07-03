function apply_plot_theme(ax, theme, varargin)
%APPLY_PLOT_THEME Apply Core v1 axis styling (non-scalogram figures).
    if nargin < 2 || isempty(theme)
        theme = plot_config('development');
    end
    if nargin < 1 || isempty(ax)
        ax = gca;
    end

    p = inputParser;
    addParameter(p, 'isScalogram', false, @islogical);
    parse(p, varargin{:});
    isScalogram = p.Results.isScalogram;

    set(ax, 'Box', 'off', 'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');

    if isScalogram
        colormap(ax, theme.scalogramColormap);
        set(ax, 'FontName', 'Times New Roman');
        return;
    end

    if strcmp(theme.mode, 'publication')
        set(ax, 'FontName', theme.axesFontName, 'FontSize', theme.axesFontSize, ...
            'FontWeight', theme.axesFontWeight, 'LineWidth', theme.axesLineWidth);
    end
end
