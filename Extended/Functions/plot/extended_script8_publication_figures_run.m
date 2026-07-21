function out = extended_script8_publication_figures_run(cohortRoot, cfg)
%EXTENDED_SCRIPT8_PUBLICATION_FIGURES_RUN Curated publication composites (Scripts 1-7 inputs).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    cfg = extended_apply_plot_cfg(cfg);
    theme = script8_theme_(cfg);
    pal = theme.palette;

    if nargin < 1 || isempty(cohortRoot)
        cohortRoot = uigetdir(pwd, 'Select cohort results folder (e.g. C57_LP)');
        if isequal(cohortRoot, 0)
            error('extended_script8_publication_figures_run:NoRoot', 'No cohort folder selected.');
        end
    end

    paths = extended_script8_resolve_paths(cohortRoot);
    if ~isempty(paths.missing)
        error('extended_script8_publication_figures_run:MissingInputs', ...
            'Missing required inputs: %s', strjoin(paths.missing, ', '));
    end

    outDirs = script8_make_output_dirs_(cohortRoot, cfg.plotMode);
    logPath = fullfile(outDirs.logs, sprintf('Script8_PublicationFigures_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() script8_fclose_(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 8: Publication figures ===\n');
    fprintf('Cohort:  %s\n', paths.cohortRoot);
    fprintf('Output:  %s\n', outDirs.root);
    fprintf('Mode:    %s | dpi=%g | ext=%s\n', cfg.plotMode, theme.dpi, theme.ext);
    script8_log_(LOG, 'Cohort: %s', paths.cohortRoot);

    data = script8_load_data_(paths, cfg);
    manifest = script8_manifest_init_();
    compositeWide = strings(0, 1);

    %% Fig 0 — Pipeline
    [manifest, wide0, stand0] = script8_build_fig00_(paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide0); %#ok<AGROW>
    script8_log_(LOG, 'Fig 0 complete');

    %% Fig 1 — Scalograms
    [manifest, wide1, tall1, stand1] = script8_build_fig01_(paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide1); %#ok<AGROW>
    script8_export_tall_copy_(tall1, wide1, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 1 complete (%d standalone)', numel(stand1));

    %% Fig 2 — Gradient
    [manifest, wide2, tall2, stand2] = script8_build_fig02_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide2); %#ok<AGROW>
    script8_export_tall_copy_(tall2, wide2, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 2 complete');

    %% Fig 3 — Transitions
    [manifest, wide3, tall3, stand3] = script8_build_fig03_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide3); %#ok<AGROW>
    script8_export_tall_copy_(tall3, wide3, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 3 complete');

    %% Fig 4 — Mechanistic
    [manifest, wide4, tall4, stand4] = script8_build_fig04_(data, paths, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide4); %#ok<AGROW>
    script8_export_tall_copy_(tall4, wide4, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 4 complete');

    %% Supplementary standalone
    [manifest, standS] = script8_build_supplementary_(data, paths, outDirs, theme, manifest);
    script8_log_(LOG, 'Supplementary complete');

    script8_write_manifest_(manifest, outDirs.manifest);
    script8_save_storyboard_(compositeWide, outDirs.storyboard);

    out = struct();
    out.cohortRoot = paths.cohortRoot;
    out.outRoot = outDirs.root;
    out.manifest = outDirs.manifest;
    out.storyboard = outDirs.storyboard;
    out.nManifestRows = height(manifest);
    out.logPath = logPath;

    fprintf('\nExtended Script 8 complete.\n');
    fprintf('  Composites: %s\n', outDirs.compositeWide);
    fprintf('  Standalone: %s\n', outDirs.standalone);
    fprintf('  Manifest:   %s\n', outDirs.manifest);
end

%% ========================================================================
function data = script8_load_data_(paths, cfg)
    data = struct();
    data.pairSummary = script8_read_sheet_(paths.lmeDescriptiveXlsx, 'CR_UR_Pairs_Summary');
    data.pairBySex = script8_read_sheet_(paths.lmeDescriptiveXlsx, 'CR_UR_Pairs_Summary_BySex');
    data.absSummary = script8_read_sheet_(paths.lmeDescriptiveXlsx, 'AbsolutePower_Summary');
    data.sexBalance = script8_read_sheet_(paths.lmeDescriptiveXlsx, 'Sex_Balance');
    data.lmeFdr = script8_read_sheet_(paths.lmeInferenceXlsx, 'LME_Inference_BH_FDR');
    data.binnedCoherence = script8_read_sheet_(paths.resyncXlsx, 'BinnedCoherence');
    data.ridgePowerFdr = script8_read_sheet_(paths.resyncXlsx, 'RidgePowerStats_BH_FDR');
    data.resyncGradient = script8_read_sheet_(paths.resyncGradientXlsx, 'Summary');
    data.phase24 = script8_read_sheet_(paths.profilesXlsx, 'PhaseCoherence_24h');
    data.ridge24 = script8_read_sheet_(paths.profilesXlsx, 'RidgePower_24h');
    data.clusterSummary = script8_read_sheet_(paths.profilesXlsx, 'ClusterSummary');
    data.retention = script8_read_sheet_(paths.gateXlsx, 'Retention_ByBand');
    data.primaryUR = string(cfg.bands.primaryUR);
    data.ppOrder = script8_pp_order_(data.pairSummary);
    data.nMice = script8_total_mice_(data.sexBalance);
end

%% Fig 0
function [manifest, widePath, standPaths] = script8_build_fig00_(paths, outDirs, theme, manifest)
    pal = theme.palette;
    standPaths = {};
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 1400 420], 'Visible', 'off');
    ax = axes(fig, 'Position', [0.05 0.22 0.9 0.68]); hold(ax, 'on');
    axis(ax, [0 10 0 1]); axis(ax, 'off');

    boxes = {
        0.3, 0.55, 'RAW activity\n(Script 2 wavelet)';
        2.0, 0.55, 'CarryForward gate\n±15% SEL_P360 (Script 4)';
        4.0, 0.55, 'Transition resync\n(Script 5)';
        6.2, 0.55, 'Gradient LME\n(Script 6)';
        4.0, 0.15, 'Phase / 24h profiles\n(Script 7)';
        2.0, 0.15, 'HSub validation\n(Core Script 1–3)';
        };
    for i = 1:size(boxes, 1)
        x = boxes{i, 1}; y = boxes{i, 2}; txt = boxes{i, 3};
        w = 1.5; h = 0.28;
        rectangle(ax, 'Position', [x - w/2, y - h/2, w, h], 'Curvature', 0.08, ...
            'FaceColor', [1 1 1], 'EdgeColor', pal.base(1, :), 'LineWidth', 1.5);
        text(ax, x, y, txt, 'HorizontalAlignment', 'center', 'FontName', theme.fontName, ...
            'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
    end
    annotation(fig, 'textbox', [0.05 0.02 0.9 0.12], 'String', ...
        'Raw UR candidates carried forward only when matched to HSub residual (±15% period tolerance).', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontName', theme.fontName, ...
        'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
    script8_arrow_(ax, 1.05, 0.55, 1.25, 0.55, pal.base(1, :));
    script8_arrow_(ax, 2.75, 0.55, 3.25, 0.55, pal.base(1, :));
    script8_arrow_(ax, 4.75, 0.55, 5.45, 0.55, pal.base(1, :));
    script8_arrow_(ax, 3.25, 0.48, 2.75, 0.28, pal.base(7, :));
    script8_arrow_(ax, 4.75, 0.42, 4.0, 0.28, pal.base(3, :));
    title(ax, 'Extended analysis pipeline (C57\_LP)', 'FontName', theme.fontName, 'FontWeight', 'bold');

    base = fullfile(outDirs.standalone, 'Fig00_Pipeline');
    standPaths = script8_export_figure_(fig, base, theme, {theme.ext, '.pdf'});
    widePath = fullfile(outDirs.compositeWide, ['Fig00_Pipeline' theme.ext]);
    copyfile(standPaths{1}, widePath, 'f');
    close(fig);
    manifest = script8_manifest_add_(manifest, 'Fig00', 'All', widePath, '16:9', ...
        'Analysis pipeline from RAW activity through CarryForward validation to transition and gradient inference.', ...
        'Script8', 'Schematic', 'Includes HSub validation branch.');
end

%% Fig 1
function [manifest, widePath, tallPath, standPaths] = script8_build_fig01_(paths, outDirs, theme, manifest)
    standDir = fullfile(outDirs.standalone, 'Fig01_Scalograms');
    extended_period_gate_ensure_dir(standDir);
    imgs = {
        paths.scalogramRawF, 'A', 'RAW stitched — Female';
        paths.scalogramRawM, 'B', 'RAW stitched — Male';
        paths.scalogramHSubF, 'C', 'HSub residual — Female';
        paths.scalogramHSubM, 'D', 'HSub residual — Male';
        };
    standPaths = cell(4, 1);
    for i = 1:4
        if i >= 3
            img = script8_crop_hsub_residual_(imgs{i, 1});
        else
            img = imread(imgs{i, 1});
        end
        fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'pixels', 'Position', [50 50 max(size(img, 2), 400) max(size(img, 1), 300)]);
        ax = axes(fig); imshow(img, 'Parent', ax); axis(ax, 'off');
        title(ax, imgs{i, 3}, 'FontName', theme.fontName, 'FontWeight', 'bold', 'Interpreter', 'none');
        script8_panel_label_(ax, imgs{i, 2}, theme);
        base = fullfile(standDir, sprintf('Fig01_%s_%s', imgs{i, 2}, strrep(imgs{i, 3}, ' ', '_')));
        standPaths{i} = script8_export_figure_(fig, base, theme, {theme.ext});
        close(fig);
        manifest = script8_manifest_add_(manifest, 'Fig01', imgs{i, 2}, standPaths{i}{1}, 'standalone', ...
            imgs{i, 3}, 'Script 3', 'Jet scalogram', '');
    end

    [widePath, tallPath] = script8_composite_grid_(standPaths, outDirs, 'Fig01_Scalograms_2x2', theme, [2 2], [16 9], [4 5]);
    manifest = script8_manifest_add_(manifest, 'Fig01', 'Composite', widePath, '16:9', ...
        'Multiscale structure across photoperiod gradient: RAW vs HSub-residual (group averages).', ...
        'Script 3', 'Jet', 'Panels A–B RAW; C–D HSub residual only.');
end

%% Fig 2
function [manifest, widePath, tallPath, standPaths] = script8_build_fig02_(data, outDirs, theme, cfg, manifest)
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig02_Gradient');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(4, 1);
    pp = data.ppOrder;
    ppLabels = arrayfun(@(x) char(script8_pp_label_(x)), pp, 'UniformOutput', false);

    %% Panel A — Delta UR-CR
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA); hold(axA, 'on');
    for b = 1:numel(data.primaryUR)
        bn = data.primaryUR(b);
        sub = data.pairSummary(string(data.pairSummary.UR_Band) == bn & string(data.pairSummary.Phase) == "All", :);
        sub = sortrows(sub, 'Photoperiod_h');
        col = script8_band_colour_(pal, bn);
        y = sub.Mean_Delta_log10;
        x = 1:height(sub);
        err = sub.SD_Delta_log10 ./ sqrt(max(data.nMice, 1));
        errorbar(axA, x, y, err, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, 'CapSize', 8);
        script8_direct_line_label_(axA, x(end), y(end), char(bn), col, theme);
        if script8_lme_pp_sig_(data.lmeFdr, bn, 'Delta')
            plot(axA, x(end), y(end), 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
        end
    end
    script8_yline_zero_(axA);
    set(axA, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axA, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axA, '\Delta(UR - CR) log_{10} power', 'FontWeight', 'bold');
    title(axA, 'CR–UR co-expression across photoperiod', 'FontWeight', 'bold');
    script8_style_axes_(axA, theme);
    script8_panel_label_(axA, 'A', theme);
    script8_add_n_annotation_(axA, data.nMice, [], theme);
    xl = xlim(axA); yl = ylim(axA);
    text(axA, xl(2), yl(1) + 0.08 * range(yl), 'L24: CR collapse \rightarrow UR > CR', ...
        'HorizontalAlignment', 'right', 'FontName', theme.fontName, 'FontSize', 9, 'Color', pal.l24);
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig02_A_Delta'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% Panel B — CR vs UR_1_3 power
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    cr = data.absSummary(string(data.absSummary.BandName) == "CR_20_28" & string(data.absSummary.Phase) == "All", :);
    ur = data.absSummary(string(data.absSummary.BandName) == "UR_1_3" & string(data.absSummary.Phase) == "All", :);
    cr = sortrows(cr, 'Photoperiod_h'); ur = sortrows(ur, 'Photoperiod_h');
    plot(axB, 1:height(cr), cr.Mean_Log10, '-s', 'Color', pal.cr, 'LineWidth', 2, 'MarkerFaceColor', pal.cr);
    plot(axB, 1:height(ur), ur.Mean_Log10, '-o', 'Color', script8_band_colour_(pal, 'UR_1_3'), 'LineWidth', 2);
    script8_direct_line_label_(axB, height(cr), cr.Mean_Log10(end), 'CR', pal.cr, theme);
    script8_direct_line_label_(axB, height(ur), ur.Mean_Log10(end), 'UR\_1\_3', script8_band_colour_(pal, 'UR_1_3'), theme);
    set(axB, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axB, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axB, 'Mean log_{10} band power', 'FontWeight', 'bold');
    title(axB, 'Absolute CR vs UR\_1\_3 power', 'FontWeight', 'bold');
    script8_style_axes_(axB, theme);
    script8_panel_label_(axB, 'B', theme);
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig02_B_AbsPower'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% Panel C — UR band heatmap (delta summary)
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axC = axes(figC);
    urBands = pal.allUR;
    M = nan(numel(urBands), numel(pp));
    for bi = 1:numel(urBands)
        sub = data.pairSummary(string(data.pairSummary.UR_Band) == urBands(bi) & string(data.pairSummary.Phase) == "All", :);
        for pi = 1:numel(pp)
            row = sub(sub.Photoperiod_h == pp(pi), :);
            if ~isempty(row), M(bi, pi) = row.Mean_Delta_log10(1); end
        end
    end
    imagesc(axC, M);
    colormap(axC, script8_diverging_cmap_());
    cb = colorbar(axC); cb.Label.String = '\Delta(UR-CR)';
    set(axC, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels, 'YTick', 1:numel(urBands), ...
        'YTickLabel', cellstr(urBands));
    title(axC, 'UR band summary (\Delta UR-CR)', 'FontWeight', 'bold');
    script8_panel_label_(axC, 'C', theme);
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, 'Fig02_C_Heatmap'), theme, {theme.ext, '.pdf'});
    close(figC);

    %% Panel D — Sex inset
    figD = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axD = axes(figD); hold(axD, 'on');
    if ~isempty(data.pairBySex)
        pooled = data.pairSummary(string(data.pairSummary.UR_Band) == "UR_1_3" & string(data.pairSummary.Phase) == "All", :);
        pooled = sortrows(pooled, 'Photoperiod_h');
        plot(axD, 1:height(pooled), pooled.Mean_Delta_log10, '--', 'Color', pal.pooled, 'LineWidth', 1.2);
        for sx = ["Female", "Male"]
            sub = data.pairBySex(string(data.pairBySex.UR_Band) == "UR_1_3" & string(data.pairBySex.Phase) == "All" & string(data.pairBySex.Sex) == sx, :);
            sub = sortrows(sub, 'Photoperiod_h');
            if sx == "Female", col = pal.female; else, col = pal.male; end
            plot(axD, 1:height(sub), sub.Mean_Delta_log10, '-o', 'Color', col, 'LineWidth', 2);
            script8_direct_line_label_(axD, height(sub), sub.Mean_Delta_log10(end), char(sx), col, theme);
        end
    end
    script8_yline_zero_(axD);
    set(axD, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axD, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axD, '\Delta(UR - CR) UR\_1\_3', 'FontWeight', 'bold');
    title(axD, 'Sex-stratified \Delta (secondary)', 'FontWeight', 'bold');
    script8_style_axes_(axD, theme);
    script8_panel_label_(axD, 'D', theme);
    standPaths{4} = script8_export_figure_(figD, fullfile(standDir, 'Fig02_D_Sex'), theme, {theme.ext, '.pdf'});
    close(figD);

    for k = 1:4
        manifest = script8_manifest_add_(manifest, 'Fig02', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            'Gradient CR–UR panel', 'Script 6', 'Tol contract', '');
    end
    [widePath, tallPath] = script8_composite_from_paths_(standPaths, outDirs, 'Fig02_Gradient_CRUR', theme, [2 2], [16 9], [4 5]);
    manifest = script8_manifest_add_(manifest, 'Fig02', 'Composite', widePath, '16:9', ...
        'Photoperiod gradient of CR–UR co-expression with sex inset.', 'Script 6', 'Tol', '');
end

%% Fig 3
function [manifest, widePath, tallPath, standPaths] = script8_build_fig03_(data, outDirs, theme, cfg, manifest)
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig03_Transitions');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(4, 1);
    pp = data.ppOrder;
    ppPlot = pp(pp <= 22);
    ppLabels = arrayfun(@(x) char(script8_pp_label_(x)), ppPlot, 'UniformOutput', false);

    grad = data.resyncGradient;
    if isempty(grad)
        grad = table();
    end

    %% A — Ridge power gradient
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA); hold(axA, 'on');
    script8_plot_transition_metric_(axA, grad, data.primaryUR, 'RidgePower_PostMinusPre', ppPlot, pal, theme);
    xlabel(axA, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axA, 'Ridge power post - pre', 'FontWeight', 'bold');
    title(axA, 'Transition-locked ridge power change', 'FontWeight', 'bold');
    script8_style_axes_(axA, theme);
    script8_panel_label_(axA, 'A', theme);
    script8_add_n_annotation_(axA, data.nMice, [], theme);
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig03_A_RidgePowerGrad'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — DeltaR gradient
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    script8_plot_transition_metric_(axB, grad, data.primaryUR, 'DeltaR', ppPlot, pal, theme);
    xlabel(axB, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axB, '\DeltaR (post - pre)', 'FontWeight', 'bold');
    title(axB, 'Transition phase resynchronisation (\DeltaR)', 'FontWeight', 'bold');
    script8_style_axes_(axB, theme);
    script8_panel_label_(axB, 'B', theme);
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig03_B_DeltaRGrad'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — Faceted coherence L12 | L22 | L24
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 420]);
    facets = pal.coherenceFacets;
    for fi = 1:numel(facets)
        axC = subplot(1, numel(facets), fi); hold(axC, 'on');
        script8_plot_coherence_facet_(axC, data.binnedCoherence, facets(fi), data.primaryUR(1), pal, theme);
        title(axC, char(script8_pp_label_(facets(fi))), 'FontWeight', 'bold');
        if fi == 1
            ylabel(axC, 'Phase coherence R', 'FontWeight', 'bold');
        end
        if fi == 2
            xlabel(axC, 'Time relative to transition (h)', 'FontWeight', 'bold');
        end
        if fi == 1
            script8_panel_label_(axC, 'C', theme);
        end
        set(axC, 'YLim', [0 pal.coherenceYMax], 'XLim', pal.coherenceXlim);
        script8_style_axes_(axC, theme);
    end
    sgtitle(figC, 'DL/LD phase coherence (UR\_1\_3, faceted by photoperiod)', 'FontWeight', 'bold', 'FontName', theme.fontName);
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, 'Fig03_C_CoherenceFacets'), theme, {theme.ext, '.pdf'});
    close(figC);

    %% D — Ridge power FDR heatmap
    figD = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axD = axes(figD);
    script8_plot_ridge_heatmap_(axD, data.ridgePowerFdr, pal, theme);
    title(axD, 'Ridge power change (FDR-significant direction)', 'FontWeight', 'bold');
    script8_panel_label_(axD, 'D', theme);
    standPaths{4} = script8_export_figure_(figD, fullfile(standDir, 'Fig03_D_Heatmap'), theme, {theme.ext, '.pdf'});
    close(figD);

    for k = 1:4
        manifest = script8_manifest_add_(manifest, 'Fig03', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            'Transition dynamics panel', 'Script 5', 'Tol DL/LD + diverging heatmap', '');
    end
    [widePath, tallPath] = script8_composite_from_paths_(standPaths, outDirs, 'Fig03_Transitions', theme, [2 2], [16 9], [4 5]);
    manifest = script8_manifest_add_(manifest, 'Fig03', 'Composite', widePath, '16:9', ...
        'Transition-locked ultradian dynamics across photoperiod.', 'Script 5', 'Tol', 'Panel C faceted L12|L22|L24.');
end

function script8_plot_transition_metric_(ax, grad, bands, metricName, ppPlot, pal, theme)
    if isempty(grad) || ~ismember('Metric', grad.Properties.VariableNames)
        return;
    end
    x = 1:numel(ppPlot);
    set(ax, 'XTick', x, 'XTickLabel', arrayfun(@(v) char(script8_pp_label_(v)), ppPlot, 'UniformOutput', false));
    script8_yline_zero_(ax);
    hold(ax, 'on');
    for tr = ["DL", "LD"]
        col = script8_transition_colour_(pal, tr);
        vals = nan(numel(ppPlot), 1);
        for pi = 1:numel(ppPlot)
            idx = grad.Photoperiod_h == ppPlot(pi) & string(grad.TransitionType) == tr & ...
                string(grad.Metric) == metricName & string(grad.BandName) == bands(1);
            if any(idx)
                vals(pi) = grad.MeanEffect(find(idx, 1));
            end
        end
        plot(ax, x, vals, '-o', 'Color', col, 'LineWidth', 2.2, 'MarkerSize', 7, 'MarkerFaceColor', col);
        if any(isfinite(vals))
            ip = find(isfinite(vals), 1, 'last');
            script8_direct_line_label_(ax, ip, vals(ip), char(tr), col, theme);
        end
    end
end

function script8_plot_coherence_facet_(ax, Bin, photoperiod_h, bandName, pal, theme)
    if isempty(Bin), return; end
    B = Bin(Bin.Photoperiod_h == photoperiod_h & string(Bin.BandName) == bandName, :);
    hold(ax, 'on');
    ylSet = false;
    for tr = ["DL", "LD"]
        Bt = B(string(B.TransitionType) == tr, :);
        if isempty(Bt), continue; end
        Bt = sortrows(Bt, 'RelBinCenter_h');
        col = script8_transition_colour_(pal, tr);
        x = Bt.RelBinCenter_h;
        y = Bt.R;
        plot(ax, x, y, '-o', 'Color', col, 'LineWidth', 2.2, 'MarkerSize', 5, 'MarkerFaceColor', col);
        if ismember('N_PhaseObs', Bt.Properties.VariableNames)
            se = 1.96 ./ sqrt(max(Bt.N_PhaseObs, 1));
            fill(ax, [x; flipud(x)], [y - se; flipud(y + se)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        end
        script8_direct_line_label_(ax, x(end), y(end), char(tr), col, theme);
        if ~ylSet
            ylim(ax, [0 pal.coherenceYMax]);
            ylSet = true;
        end
        script8_shade_transition_(ax, pal.coherenceXlim, tr == "DL", pal);
    end
    xline(ax, 0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    script8_style_axes_(ax, theme);
end

function script8_plot_ridge_heatmap_(ax, RP, pal, theme)
    if isempty(RP), return; end
    bands = pal.allUR;
    pp = unique(double(RP.Photoperiod_h)); pp = sort(pp(pp <= 22));
    M = nan(numel(bands), numel(pp));
    for bi = 1:numel(bands)
        for pi = 1:numel(pp)
            idx = RP.Photoperiod_h == pp(pi) & string(RP.BandName) == bands(bi) & string(RP.TransitionType) == "LD";
            if any(idx) && ismember('MeanDifference_PostMinusPre', RP.Properties.VariableNames)
                M(bi, pi) = RP.MeanDifference_PostMinusPre(find(idx, 1));
            end
        end
    end
    imagesc(ax, M);
    colormap(ax, script8_diverging_cmap_());
    colorbar(ax);
    set(ax, 'XTick', 1:numel(pp), 'XTickLabel', arrayfun(@(x) char(script8_pp_label_(x)), pp, 'UniformOutput', false), ...
        'YTick', 1:numel(bands), 'YTickLabel', cellstr(bands));
    script8_style_axes_(ax, theme);
end

%% Fig 4
function [manifest, widePath, tallPath, standPaths] = script8_build_fig04_(data, paths, outDirs, theme, cfg, manifest)
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig04_Mechanistic');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(4, 1);

    %% A — Cluster histogram
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA); hold(axA, 'on');
    CS = data.clusterSummary;
    if ~isempty(CS) && ismember('PeriodCentre_h', CS.Properties.VariableNames)
        histogram(axA, CS.PeriodCentre_h, 'BinWidth', 0.25, 'FaceColor', pal.base(2, :), 'EdgeColor', 'k');
        xline(axA, 2, '--', 'Color', script8_band_colour_(pal, 'UR_1_3'), 'LineWidth', 1.2);
        xline(axA, 4.5, '--', 'Color', script8_band_colour_(pal, 'UR_3_6'), 'LineWidth', 1.2);
    end
    xlabel(axA, 'Validated ridge period (h)', 'FontWeight', 'bold');
    ylabel(axA, 'Count', 'FontWeight', 'bold');
    title(axA, 'Validated period clusters', 'FontWeight', 'bold');
    script8_style_axes_(axA, theme);
    script8_panel_label_(axA, 'A', theme);
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig04_A_Clusters'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — 24h phase L12 vs L24
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    script8_plot_24h_profile_(axB, data.phase24, 'R', [12 24], pal, theme, 'Phase coherence R');
    title(axB, '24h phase coherence (UR\_1\_3)', 'FontWeight', 'bold');
    script8_panel_label_(axB, 'B', theme);
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig04_B_Phase24h'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — 24h ridge power
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axC = axes(figC); hold(axC, 'on');
    script8_plot_24h_profile_(axC, data.ridge24, 'MeanRidgePower_log10', [12 24], pal, theme, 'Ridge power (log_{10})');
    title(axC, '24h ridge power (UR\_1\_3)', 'FontWeight', 'bold');
    script8_panel_label_(axC, 'C', theme);
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, 'Fig04_C_Ridge24h'), theme, {theme.ext, '.pdf'});
    close(figC);

    %% D — Retention by band
    figD = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axD = axes(figD); hold(axD, 'on');
    R = data.retention;
    if ~isempty(R)
        bn = string(R.BandName);
        retCol = R.Properties.VariableNames{contains(R.Properties.VariableNames, 'RetentionPct_Elig', 'IgnoreCase', true)};
        y = double(R.(retCol));
        cols = arrayfun(@(i) script8_band_colour_(pal, bn(i)), 1:numel(bn), 'UniformOutput', false);
        for i = 1:numel(y)
            bar(axD, i, y(i), 0.7, 'FaceColor', cols{i}, 'EdgeColor', 'k');
        end
        set(axD, 'XTick', 1:numel(bn), 'XTickLabel', cellstr(bn));
        ylim(axD, [0 105]);
    end
    ylabel(axD, 'Validated retention (% eligible Raw)', 'FontWeight', 'bold');
    title(axD, 'CarryForward retention by band', 'FontWeight', 'bold');
    script8_style_axes_(axD, theme);
    script8_panel_label_(axD, 'D', theme);
    standPaths{4} = script8_export_figure_(figD, fullfile(standDir, 'Fig04_D_Retention'), theme, {theme.ext, '.pdf'});
    close(figD);

    for k = 1:4
        manifest = script8_manifest_add_(manifest, 'Fig04', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            'Mechanistic / validation panel', 'Scripts 4/7', 'Tol', '');
    end
    [widePath, tallPath] = script8_composite_from_paths_(standPaths, outDirs, 'Fig04_Mechanistic', theme, [2 2], [16 9], [4 5]);
    manifest = script8_manifest_add_(manifest, 'Fig04', 'Composite', widePath, '16:9', ...
        'Mechanistic depth: clusters, 24h profiles, CarryForward retention.', 'Scripts 4/7', 'Tol', '');
end

function script8_plot_24h_profile_(ax, T, yCol, photos, pal, theme, yLabel)
    if isempty(T) || ~ismember(yCol, T.Properties.VariableNames)
        return;
    end
    if ismember('BandName', T.Properties.VariableNames)
        T = T(string(T.BandName) == "UR_1_3", :);
    end
    hold(ax, 'on');
    cols = {pal.l12, pal.l24};
    labs = {script8_pp_label_(photos(1)), script8_pp_label_(photos(2))};
    for i = 1:numel(photos)
        sub = T(T.Photoperiod_h == photos(i), :);
        if isempty(sub), continue; end
        xCol = 'ZTBinCenter_h';
        if ~ismember(xCol, sub.Properties.VariableNames)
            xCol = 'RelBinCenter_h';
        end
        if ~ismember(xCol, sub.Properties.VariableNames), continue; end
        sub = sortrows(sub, xCol);
        subG = groupsummary(sub, xCol, 'mean', yCol);
        meanCol = "mean_" + string(yCol);
        if ~ismember(meanCol, subG.Properties.VariableNames)
            vn = subG.Properties.VariableNames(startsWith(subG.Properties.VariableNames, 'mean_'));
            if ~isempty(vn), meanCol = string(vn{1}); end
        end
        plot(ax, subG.(xCol), subG.(meanCol), '-', 'Color', cols{i}, 'LineWidth', 2.2);
        script8_direct_line_label_(ax, subG.(xCol)(end), subG.(meanCol)(end), char(labs{i}), cols{i}, theme);
    end
    xlabel(ax, 'ZT (h)', 'FontWeight', 'bold');
    ylabel(ax, yLabel, 'FontWeight', 'bold');
    script8_style_axes_(ax, theme);
end

%% Supplementary
function [manifest, standPaths] = script8_build_supplementary_(data, paths, outDirs, theme, manifest)
    standDir = fullfile(outDirs.standalone, 'Supplementary');
    extended_period_gate_ensure_dir(standDir);
    standPaths = {};
    coexpPath = fullfile(paths.script3Root, 'Figures', 'Coexpression_CR_UR_RAW.jpeg');
    if isfile(coexpPath)
        dest = fullfile(standDir, ['Supp_Coexpression_RAW' theme.ext]);
        copyfile(coexpPath, dest, 'f');
        standPaths{end + 1} = dest; %#ok<AGROW>
        manifest = script8_manifest_add_(manifest, 'Supp', 'Coexpression', dest, 'standalone', ...
            'Core exploratory RAW co-expression (Script 3).', 'Script 3', 'Collaborator', '');
    end
end

%% Composite helpers
function [widePath, tallPath] = script8_composite_from_paths_(standPaths, outDirs, name, theme, gridShape, wideAR, tallAR) %#ok<INUSD>
    imgs = cell(numel(standPaths), 1);
    for i = 1:numel(standPaths)
        p = standPaths{i};
        if iscell(p), p = p{1}; end
        if isfile(p), imgs{i} = imread(p); end
    end
    [widePath, tallPath] = script8_composite_grid_(imgs, outDirs, name, theme, gridShape, wideAR, tallAR);
end

function [widePath, tallPath] = script8_composite_grid_(items, outDirs, name, theme, gridShape, wideAR, tallAR) %#ok<INUSD>
    n = gridShape(1) * gridShape(2);
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 900]);
    tl = tiledlayout(figW, gridShape(1), gridShape(2), 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:min(n, numel(items))
        nexttile(tl, i);
        img = items{i};
        if iscell(img), img = imread(img{1}); end
        if ~isempty(img), imshow(img); axis off; end
    end
    widePath = fullfile(outDirs.compositeWide, [name theme.ext]);
    exportgraphics(figW, widePath, 'Resolution', theme.dpi);
    close(figW);

    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1250]);
    tl2 = tiledlayout(figT, gridShape(1), gridShape(2), 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:min(n, numel(items))
        nexttile(tl2, i);
        img = items{i};
        if iscell(img), img = imread(img{1}); end
        if ~isempty(img), imshow(img); axis off; end
    end
    tallPath = fullfile(outDirs.compositeTall, [name theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
end

function script8_export_tall_copy_(tallPath, widePath, tallDir)
    if isfile(tallPath)
        return;
    end
    if isfile(widePath)
        copyfile(widePath, fullfile(tallDir, ['TALL_' basename_(widePath)]), 'f');
    end
end

function b = basename_(p)
    [~, b, e] = fileparts(p);
    b = [b e];
end

%% Theme / IO / manifest helpers
function theme = script8_theme_(cfg)
    pal = extended_tol_bright_palette();
    theme = struct('palette', pal, 'fontName', pal.fontName, ...
        'dpi', cfg.plot.saveDpi, 'ext', cfg.plot.figExt, 'plotMode', string(cfg.plotMode));
end

function outDirs = script8_make_output_dirs_(cohortRoot, plotMode)
    modeLabel = "Development";
    if strcmpi(plotMode, 'publication'), modeLabel = "Publication"; end
    outRoot = fullfile(cohortRoot, ['Script8_PublicationFigures_' char(modeLabel)]);
    outDirs = struct('root', outRoot, ...
        'compositeWide', fullfile(outRoot, 'CompositePanels', 'Wide_16x9'), ...
        'compositeTall', fullfile(outRoot, 'CompositePanels', 'Tall_4x5'), ...
        'standalone', fullfile(outRoot, 'Standalone'), ...
        'logs', fullfile(outRoot, 'Logs'), ...
        'manifest', fullfile(outRoot, 'Manifest.xlsx'), ...
        'storyboard', fullfile(outRoot, 'Storyboard_AllFigures.pdf'));
    extended_period_gate_ensure_dir(outDirs.root);
    extended_period_gate_ensure_dir(outDirs.compositeWide);
    extended_period_gate_ensure_dir(outDirs.compositeTall);
    extended_period_gate_ensure_dir(outDirs.standalone);
    extended_period_gate_ensure_dir(outDirs.logs);
end

function manifest = script8_manifest_init_()
    manifest = table('Size', [0 8], ...
        'VariableTypes', repmat({'string'}, 1, 8), ...
        'VariableNames', {'Figure', 'Panel', 'File', 'Aspect', 'Caption', 'SourceScript', 'ColourKey', 'Notes'});
end

function manifest = script8_manifest_add_(manifest, figId, panelId, filePath, aspect, caption, sourceScript, colourKey, notes)
    if nargin < 9, notes = ""; end
    manifest = [manifest; {string(figId), string(panelId), string(filePath), string(aspect), ...
        string(caption), string(sourceScript), string(colourKey), string(notes)}]; %#ok<AGROW>
end

function script8_write_manifest_(manifest, xlsxPath)
    writetable(manifest, xlsxPath, 'Sheet', 'Manifest');
    writetable(script8_colour_key_table_(), xlsxPath, 'Sheet', 'ColourKey');
end

function T = script8_colour_key_table_()
    roles = {
        'DL (lights-on)', '#4477AA', 'Transition panels only';
        'LD (lights-off)', '#EE6677', 'Transition panels only';
        'UR_1_3', '#66CCEE', 'Band panels';
        'UR_3_6', '#AA3377', 'Band panels';
        'CR_20_28', '#BBBBBB', 'CR traces / dark schematic';
        'Female', '#228833', 'Sex inset only';
        'Male', '#CCBB44', 'Sex inset only';
        'L12 (24h compare)', '#4477AA', 'Fig 4B/C only';
        'L24 (24h compare)', '#AA3377', 'Fig 4B/C only';
        'Scalograms', 'Jet', 'Fig 1 only'};
    T = cell2table(roles, 'VariableNames', {'Role', 'Colour', 'Usage'});
end

function script8_save_storyboard_(compositeFiles, pdfPath)
    compositeFiles = compositeFiles(isfile(compositeFiles));
    if isempty(compositeFiles), return; end
    for i = 1:numel(compositeFiles)
        img = imread(compositeFiles{i});
        fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 size(img, 2) size(img, 1)]);
        ax = axes(fig, 'Position', [0 0 1 1]); %#ok<LAXES>
        imshow(img, 'Parent', ax); axis(ax, 'off');
        if i == 1
            exportgraphics(fig, pdfPath, 'ContentType', 'image');
        else
            exportgraphics(fig, pdfPath, 'ContentType', 'image', 'Append', true);
        end
        close(fig);
    end
end

function outPaths = script8_export_figure_(fig, basePath, theme, formats)
    outPaths = cell(numel(formats), 1);
    extended_period_gate_ensure_dir(fileparts(basePath));
    for i = 1:numel(formats)
        ext = formats{i};
        if ~startsWith(ext, '.'), ext = ['.' ext]; end
        outFile = [basePath ext];
        exportgraphics(fig, outFile, 'Resolution', theme.dpi);
        outPaths{i} = outFile;
    end
end

function T = script8_read_sheet_(path, sheet)
    try
        T = readtable(path, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    catch
        T = table();
    end
end

function ppOrder = script8_pp_order_(T)
    if isempty(T) || ~ismember('Photoperiod_h', T.Properties.VariableNames)
        ppOrder = [12 14 16 18 20 22 24];
        return;
    end
    ppOrder = unique(double(T.Photoperiod_h));
    ppOrder = sort(ppOrder(isfinite(ppOrder)));
end

function n = script8_total_mice_(sexBalance)
    n = 24;
    if isempty(sexBalance) || ~ismember('N_Mice', sexBalance.Properties.VariableNames)
        return;
    end
    if ismember('Sex', sexBalance.Properties.VariableNames)
        G = findgroups(sexBalance.Sex);
        mx = splitapply(@max, sexBalance.N_Mice, G);
        n = sum(mx);
    else
        n = max(sexBalance.N_Mice);
    end
end

function lbl = script8_pp_label_(pp)
    lbl = "L" + string(round(double(pp)));
end

function rgb = script8_band_colour_(pal, bandName)
    key = char(string(bandName));
    if isKey(pal.band, key), rgb = pal.band(key); else, rgb = pal.base(7, :); end
end

function rgb = script8_transition_colour_(pal, transitionType)
    t = upper(string(transitionType));
    if t == "DL", rgb = pal.dl; elseif t == "LD", rgb = pal.ld; else, rgb = pal.base(7, :); end
end

function tf = script8_lme_pp_sig_(lmeTable, bandName, metricClass)
    tf = false;
    if isempty(lmeTable), return; end
    sub = lmeTable(string(lmeTable.BandName) == bandName & ...
        string(lmeTable.MetricClass) == metricClass & string(lmeTable.Term) == "Photoperiod_h", :);
    if ~isempty(sub) && ismember('Significant_BH', sub.Properties.VariableNames)
        tf = any(sub.Significant_BH);
    end
end

function script8_style_axes_(ax, theme)
    set(ax, 'Box', 'off', 'TickDir', theme.palette.tickDir, 'FontName', theme.fontName, 'LineWidth', 1);
end

function script8_panel_label_(ax, labelChar, theme)
    text(ax, 0.02, 0.98, labelChar, 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 14, 'FontName', theme.fontName, 'VerticalAlignment', 'top');
end

function script8_direct_line_label_(ax, x, y, txt, colour, theme)
    text(ax, x, y, txt, 'Color', colour, 'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 10);
end

function script8_add_n_annotation_(ax, nMice, ~, theme)
    text(ax, 0.98, 0.02, sprintf('n=%d mice', nMice), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontName', theme.fontName, 'FontSize', 9);
end

function script8_yline_zero_(ax)
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
end

function script8_shade_transition_(ax, xlimRange, isDL, pal)
    yl = ylim(ax);
    if isDL
        patch(ax, [xlimRange(1) 0 0 xlimRange(1)], [yl(1) yl(1) yl(2) yl(2)], pal.cr, ...
            'FaceAlpha', 0.07, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    else
        patch(ax, [0 xlimRange(2) xlimRange(2) 0], [yl(1) yl(1) yl(2) yl(2)], pal.cr, ...
            'FaceAlpha', 0.07, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function img = script8_crop_hsub_residual_(imgPath)
    img = imread(imgPath);
    img = img(:, round(size(img, 2) / 2) + 1:end, :);
end

function script8_arrow_(ax, x1, y1, x2, y2, col)
    annotation(ax.Parent, 'arrow', ...
        script8_norm_ann_(ax, [x1 x2]), script8_norm_ann_y_(ax, [y1 y2]), ...
        'Color', col, 'LineWidth', 1.5, 'HeadWidth', 8);
end

function xn = script8_norm_ann_(ax, x)
    xl = xlim(ax); xn = (x - xl(1)) / range(xl) * 0.9 + 0.05;
end

function yn = script8_norm_ann_y_(ax, y)
    yl = ylim(ax); yn = (y - yl(1)) / range(yl) * 0.9 + 0.05;
end

function cmap = script8_diverging_cmap_()
    pal = extended_tol_bright_palette();
    n = 32;
    low = [linspace(pal.base(1, 1), 1, n)', linspace(pal.base(1, 2), 1, n)', linspace(pal.base(1, 3), 1, n)'];
    hi = [linspace(1, pal.base(5, 1), n)', linspace(1, pal.base(5, 2), n)', linspace(1, pal.base(5, 3), n)'];
    cmap = [low; hi];
end

function script8_log_(LOG, fmt, varargin)
    if ~isempty(LOG) && LOG > 0
        fprintf(LOG, ['[' datestr(now, 'yyyy-mm-dd HH:MM:SS') '] ' fmt '\n'], varargin{:});
    end
end

function script8_fclose_(LOG)
    if ~isempty(LOG) && LOG > 0
        fclose(LOG);
    end
end
