function theme = plot_theme_ensure_scalogram(theme)
%PLOT_THEME_ENSURE_SCALOGRAM Ensure theme.scalogram exists (callable from any script).
    if nargin < 1 || isempty(theme)
        theme = plot_config('development');
        return;
    end
    if isfield(theme, 'scalogram') && isstruct(theme.scalogram) && isfield(theme.scalogram, 'format')
        return;
    end
    cfg = core_defaults();
    mode = 'development';
    if isfield(theme, 'mode') && ~isempty(theme.mode)
        mode = lower(char(theme.mode));
    end
    if isfield(cfg.plot.scalogram, mode)
        theme.scalogram = cfg.plot.scalogram.(mode);
    else
        theme.scalogram = cfg.plot.scalogram.development;
    end
    if ~isfield(theme, 'scalogramColormap')
        theme.scalogramColormap = theme.scalogram.colormap;
    end
end
