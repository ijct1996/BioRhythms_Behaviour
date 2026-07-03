function export_figure(fig, outPath, theme, varargin)
%EXPORT_FIGURE Save figure with mode-appropriate DPI/format.
%
%   export_figure(gcf, outFile, plot_config('development'))
%   export_figure(gcf, outFile, plot_config('publication'), 'isScalogram', true)

    if nargin < 3 || isempty(theme)
        theme = plot_config('development');
    end
    theme = plot_theme_ensure_scalogram(theme);

    p = inputParser;
    addParameter(p, 'isScalogram', false, @islogical);
    addParameter(p, 'dpiOverride', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'formatOverride', '', @ischar);
    parse(p, varargin{:});

    [outDir, ~, ~] = fileparts(outPath);
    if ~isempty(outDir)
        ensure_dir(outDir);
    end

    if p.Results.isScalogram
        dpi = theme.scalogram.dpi;
        fmt = lower(theme.scalogram.format);
    else
        dpi = theme.dpi;
        fmt = lower(theme.exportFormat);
    end
    if ~isempty(p.Results.dpiOverride)
        dpi = p.Results.dpiOverride;
    end
    if ~isempty(p.Results.formatOverride)
        fmt = lower(p.Results.formatOverride);
    end

    switch fmt
        case 'jpeg'
            print(fig, outPath, '-djpeg', sprintf('-r%d', dpi));
        case 'png'
            print(fig, outPath, '-dpng', sprintf('-r%d', dpi));
        otherwise
            saveas(fig, outPath);
    end
end
