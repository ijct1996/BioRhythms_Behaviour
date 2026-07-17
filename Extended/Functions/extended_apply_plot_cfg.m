function cfg = extended_apply_plot_cfg(cfg)
%EXTENDED_APPLY_PLOT_CFG Set dpi/format from cfg.plotMode (development | publication).

    if nargin < 1 || isempty(cfg)
        cfg = extended_defaults();
        return;
    end

    if ~isfield(cfg, 'plotMode') || isempty(cfg.plotMode)
        cfg.plotMode = 'development';
    end

    if ~isfield(cfg, 'plot'), cfg.plot = struct(); end

    coreCfg = core_defaults();
    switch lower(string(cfg.plotMode))
        case "publication"
            cfg.plot.saveDpi = coreCfg.plot.publication.dpi;
            cfg.plot.figExt = ['.' coreCfg.plot.publication.format];
            cfg.plot.saveTiff = true;
        otherwise
            cfg.plot.saveDpi = coreCfg.plot.development.dpi;
            cfg.plot.figExt = ['.' coreCfg.plot.development.format];
            cfg.plot.saveTiff = false;
    end
end
