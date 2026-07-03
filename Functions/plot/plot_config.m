function theme = plot_config(mode)
%PLOT_CONFIG Development vs publication theme for Core v1.
%
%   theme = plot_config('development')   % default
%   theme = plot_config('publication')

    cfg = core_defaults();
    pal = collaborator_palette();

    if nargin < 1 || isempty(mode)
        mode = 'development';
    end

    theme.mode = lower(mode);
    theme.scalogramColormap = 'jet';

    switch theme.mode
        case 'publication'
            theme.fontName = pal.publication.fontName;
            theme.dpi = cfg.plot.publication.dpi;
            theme.exportFormat = cfg.plot.publication.format;
            theme.useCollaboratorPalette = true;
            theme.axesFontName = pal.axes.fontName;
            theme.axesFontSize = pal.axes.fontSize;
            theme.axesFontWeight = pal.axes.fontWeight;
            theme.axesLineWidth = pal.axes.lineWidth;
        otherwise
            theme.fontName = 'default';
            theme.dpi = cfg.plot.development.dpi;
            theme.exportFormat = cfg.plot.development.format;
            theme.useCollaboratorPalette = false;
            theme.axesFontName = 'default';
            theme.axesFontSize = 10;
            theme.axesFontWeight = 'normal';
            theme.axesLineWidth = 0.5;
    end

    theme.scalogram = cfg.plot.scalogram.(theme.mode);
    theme.palette = pal;
    theme = plot_theme_ensure_scalogram(theme);
end
