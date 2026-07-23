function out = extended_script8_publication_figures_run(cohortRoot, cfg)
%EXTENDED_SCRIPT8_PUBLICATION_FIGURES_RUN Curated publication composites (Scripts 1-7 inputs).
%
%   Five main figures:
%     Fig01 High-level methods schematic (5 steps)
%     Fig02 RAW + HSub scalograms + CarryForward retention (A-E)
%     Fig03 Photoperiod gradient CR–UR (A-C; sex in Supplementary)
%     Fig04 Transition coherence hero (C then A|B)
%     Fig05 Cluster-anchored activity resync at LD/DL

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    cfg = extended_apply_plot_cfg(cfg);
    theme = script8_theme_(cfg);

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

    %% Fig01 — Methods overview
    [manifest, wide1, ~] = script8_build_fig01_(paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide1); %#ok<AGROW>
    script8_log_(LOG, 'Fig 01 complete');

    %% Fig02 — RAW + HSub + retention
    [manifest, wide2, tall2, stand2] = script8_build_fig02_(data, paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide2); %#ok<AGROW>
    script8_export_tall_copy_(tall2, wide2, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 02 complete (%d standalone)', numel(stand2));

    %% Fig03 — Gradient CR–UR (A–C; sex → Supplementary)
    [manifest, wide3, tall3] = script8_build_fig03_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide3); %#ok<AGROW>
    script8_export_tall_copy_(tall3, wide3, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 03 complete');

    %% Fig04 — Transitions (C hero, then A|B)
    [manifest, wide4, tall4] = script8_build_fig04_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide4); %#ok<AGROW>
    script8_export_tall_copy_(tall4, wide4, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 04 complete');

    %% Fig05 — Cluster + transition-aligned activity resync
    [manifest, wide5, tall5] = script8_build_fig05_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide5); %#ok<AGROW>
    script8_export_tall_copy_(tall5, wide5, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 05 complete');

    %% Supplementary
    [manifest, ~] = script8_build_supplementary_(data, paths, outDirs, theme, manifest);
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
    data.activityZT = script8_read_sheet_(paths.profilesXlsx, 'ActivityComponent_24h');
    data.primaryUR = string(cfg.bands.primaryUR);
    data.ppOrder = script8_pp_order_(data.pairSummary);
    data.nMice = script8_total_mice_(data.sexBalance);
end

%% Fig01 — High-level methods schematic (publication narrative)
function [manifest, widePath, standPaths] = script8_build_fig01_(paths, outDirs, theme, manifest) %#ok<INUSD>
    pal = theme.palette;
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 1400 380], 'Visible', 'off');
    ax = axes(fig, 'Position', [0.04 0.20 0.92 0.68]); hold(ax, 'on');
    axis(ax, [0 10 0 1]); axis(ax, 'off');

    steps = {
        1.0, sprintf('Raw\nactivity');
        3.0, sprintf('Wavelet');
        5.0, sprintf('Harmonic\nsubtraction');
        7.0, sprintf('CR–UR\nco-expression');
        9.0, sprintf('Ultradian\nphase coherence');
        };
    w = 1.55; h = 0.38;
    for i = 1:size(steps, 1)
        x = steps{i, 1}; txt = steps{i, 2};
        rectangle(ax, 'Position', [x - w/2, 0.45 - h/2, w, h], 'Curvature', 0.08, ...
            'FaceColor', [1 1 1], 'EdgeColor', pal.base(1, :), 'LineWidth', 1.6);
        text(ax, x, 0.45, txt, 'HorizontalAlignment', 'center', 'FontName', theme.fontName, ...
            'FontWeight', 'bold', 'FontSize', 13, 'Interpreter', 'none');
        text(ax, x, 0.72, sprintf('%d', i), 'HorizontalAlignment', 'center', ...
            'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 14, 'Color', pal.base(1, :));
        if i < size(steps, 1)
            text(ax, x + 1.0, 0.45, '\rightarrow', 'Color', pal.base(1, :), 'FontSize', 24, ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end

    text(ax, 5, 0.08, ...
        'Extended validation: Raw UR ridge candidates carried forward only when period-matched to HSub residual (+/-15% SEL_P360).', ...
        'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontSize', 11, ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    title(ax, 'Analysis overview', 'FontName', theme.fontName, 'FontWeight', 'bold');

    base = fullfile(outDirs.standalone, 'Fig01_MethodsOverview');
    standPaths = script8_export_figure_(fig, base, theme, {theme.ext, '.pdf'});
    widePath = fullfile(outDirs.compositeWide, ['Fig01_MethodsOverview' theme.ext]);
    copyfile(standPaths{1}, widePath, 'f');
    close(fig);
    manifest = script8_manifest_add_(manifest, 'Fig01', 'All', widePath, '16:9', ...
        'High-level analysis steps: raw activity, wavelet, harmonic subtraction, CR–UR co-expression, ultradian phase coherence.', ...
        'Script8', 'Schematic', 'Script numbers omitted from boxes; see caption.');
end

%% Fig02 — RAW A|B, HSub C|D, retention E
function [manifest, widePath, tallPath, standPaths] = script8_build_fig02_(data, paths, outDirs, theme, manifest)
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig02_Scalograms');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(5, 1);

    panels = {
        paths.scalogramRawF,  'A', 'RAW stitched — Female',   false;
        paths.scalogramRawM,  'B', 'RAW stitched — Male',     false;
        paths.scalogramHSubF, 'C', 'HSub residual — Female',  true;
        paths.scalogramHSubM, 'D', 'HSub residual — Male',    true;
        };
    for i = 1:4
        if ~isfile(panels{i, 1})
            error('extended_script8_publication_figures_run:MissingScalogram', ...
                'Missing scalogram image: %s', panels{i, 1});
        end
        if panels{i, 4}
            img = script8_crop_hsub_residual_(panels{i, 1});
        else
            img = imread(panels{i, 1});
        end
        fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'pixels', ...
            'Position', [50 50 max(size(img, 2), 400) max(size(img, 1), 300)]);
        ax = axes(fig); imshow(img, 'Parent', ax); axis(ax, 'off');
        title(ax, panels{i, 3}, 'FontName', theme.fontName, 'FontWeight', 'bold', 'Interpreter', 'none');
        script8_panel_label_(ax, panels{i, 2}, theme);
        base = fullfile(standDir, sprintf('Fig02_%s_%s', panels{i, 2}, ...
            regexprep(panels{i, 3}, '[^A-Za-z0-9]+', '_')));
        standPaths{i} = script8_export_figure_(fig, base, theme, {theme.ext});
        close(fig);
        note = 'RAW circadian context in main figure.';
        if panels{i, 4}, note = 'HSub residual half-crop.'; end
        manifest = script8_manifest_add_(manifest, 'Fig02', panels{i, 2}, standPaths{i}{1}, 'standalone', ...
            panels{i, 3}, 'Script 3', 'Jet scalogram', note);
    end

    %% Panel E — CarryForward retention
    figE = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 720 420]);
    axE = axes(figE); hold(axE, 'on');
    R = data.retention;
    if ~isempty(R)
        bn = string(R.BandName);
        retCol = R.Properties.VariableNames{contains(R.Properties.VariableNames, 'RetentionPct_Elig', 'IgnoreCase', true)};
        y = double(R.(retCol));
        for i = 1:numel(y)
            bar(axE, i, y(i), 0.7, 'FaceColor', script8_band_colour_(pal, bn(i)), 'EdgeColor', 'k');
        end
        set(axE, 'XTick', 1:numel(bn), 'XTickLabel', arrayfun(@(b) script8_band_display_(b, 'tex'), bn, 'UniformOutput', false));
        axE.TickLabelInterpreter = 'tex';
        ylim(axE, [0 105]);
    end
    ylabel(axE, 'Validated retention (% eligible Raw)', 'FontWeight', 'bold');
    title(axE, 'CarryForward retention by band', 'FontWeight', 'bold');
    script8_style_axes_(axE, theme);
    script8_panel_label_(axE, 'E', theme);
    standPaths{5} = script8_export_figure_(figE, fullfile(standDir, 'Fig02_E_Retention'), theme, {theme.ext, '.pdf'});
    close(figE);
    manifest = script8_manifest_add_(manifest, 'Fig02', 'E', standPaths{5}{1}, 'standalone', ...
        'CarryForward retention by band', 'Script 4', 'Tol bands', '');

    [widePath, tallPath] = script8_composite_fig02_(standPaths, outDirs, theme);
    manifest = script8_manifest_add_(manifest, 'Fig02', 'Composite', widePath, '16:9', ...
        'RAW and HSub residual scalograms with CarryForward retention QC.', ...
        'Scripts 3/4', 'Jet + Tol', 'Tall primary: A|B, C|D, E.');
end

function [widePath, tallPath] = script8_composite_fig02_(standPaths, outDirs, theme)
    imgs = cell(5, 1);
    for i = 1:5
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script8_trim_image_whitespace_(imread(p));
    end
    % Wide: same stacking as tall (A|B / C|D / E)
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 1100]);
    tl = tiledlayout(figW, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1); imshow(imgs{1}); axis off;
    nexttile(tl, 2); imshow(imgs{2}); axis off;
    nexttile(tl, 3); imshow(imgs{3}); axis off;
    nexttile(tl, 4); imshow(imgs{4}); axis off;
    nexttile(tl, 5, [1 2]); imshow(imgs{5}); axis off;
    widePath = fullfile(outDirs.compositeWide, ['Fig02_Scalograms_Retention' theme.ext]);
    exportgraphics(figW, widePath, 'Resolution', theme.dpi);
    close(figW);

    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1400]);
    tl2 = tiledlayout(figT, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl2, 1); imshow(imgs{1}); axis off;
    nexttile(tl2, 2); imshow(imgs{2}); axis off;
    nexttile(tl2, 3); imshow(imgs{3}); axis off;
    nexttile(tl2, 4); imshow(imgs{4}); axis off;
    nexttile(tl2, 5, [1 2]); imshow(imgs{5}); axis off;
    tallPath = fullfile(outDirs.compositeTall, ['Fig02_Scalograms_Retention' theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
end

%% Fig03 — Gradient CR–UR (A–C main; D sex also written for Supplementary reuse)
function [manifest, widePath, tallPath, standPaths] = script8_build_fig03_(data, outDirs, theme, cfg, manifest) %#ok<INUSD>
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig03_Gradient');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(4, 1);
    pp = data.ppOrder;
    ppLabels = arrayfun(@(x) char(script8_pp_label_(x)), pp, 'UniformOutput', false);

    %% A — UR band heatmap
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA);
    urBands = pal.allUR;
    M = nan(numel(urBands), numel(pp));
    for bi = 1:numel(urBands)
        sub = data.pairSummary(string(data.pairSummary.UR_Band) == urBands(bi) & string(data.pairSummary.Phase) == "All", :);
        for pi = 1:numel(pp)
            row = sub(sub.Photoperiod_h == pp(pi), :);
            if ~isempty(row), M(bi, pi) = row.Mean_Delta_log10(1); end
        end
    end
    imagesc(axA, M);
    colormap(axA, script8_diverging_cmap_());
    cb = colorbar(axA); cb.Label.String = '\Delta(UR-CR)';
    yLab = arrayfun(@(b) script8_band_display_(b, 'tex'), urBands, 'UniformOutput', false);
    set(axA, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels, 'YTick', 1:numel(urBands), ...
        'YTickLabel', yLab, 'TickLabelInterpreter', 'tex');
    title(axA, 'UR band summary (\Delta UR-CR)', 'FontWeight', 'bold');
    script8_panel_label_(axA, 'A', theme);
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig03_A_Heatmap'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — Absolute power UR 1–3 + UR 3–6 (+ CR reference)
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    cr = data.absSummary(string(data.absSummary.BandName) == "CR_20_28" & string(data.absSummary.Phase) == "All", :);
    cr = sortrows(cr, 'Photoperiod_h');
    crErr = cr.SD_Log10 ./ sqrt(max(data.nMice, 1));
    errorbar(axB, 1:height(cr), cr.Mean_Log10, crErr, '-s', 'Color', pal.cr, 'LineWidth', 1.5, ...
        'MarkerFaceColor', pal.cr, 'CapSize', 6);
    script8_direct_line_label_(axB, height(cr), cr.Mean_Log10(end), script8_band_display_('CR_20_28', 'tex'), pal.cr, theme, true);
    for b = 1:numel(data.primaryUR)
        bn = data.primaryUR(b);
        ur = data.absSummary(string(data.absSummary.BandName) == bn & string(data.absSummary.Phase) == "All", :);
        ur = sortrows(ur, 'Photoperiod_h');
        col = script8_band_colour_(pal, bn);
        urErr = ur.SD_Log10 ./ sqrt(max(data.nMice, 1));
        errorbar(axB, 1:height(ur), ur.Mean_Log10, urErr, '-o', 'Color', col, 'LineWidth', 2, ...
            'MarkerFaceColor', col, 'CapSize', 8);
        script8_direct_line_label_(axB, height(ur), ur.Mean_Log10(end), script8_band_display_(bn, 'tex'), col, theme, true);
    end
    set(axB, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axB, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axB, 'Mean log_{10} band power', 'FontWeight', 'bold');
    title(axB, 'Absolute power: UR 1–3 and UR 3–6 (CR reference)', 'FontWeight', 'bold');
    script8_style_axes_(axB, theme);
    script8_panel_label_(axB, 'B', theme);
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig03_B_AbsPower'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — Delta co-expression both primary URs
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axC = axes(figC); hold(axC, 'on');
    for b = 1:numel(data.primaryUR)
        bn = data.primaryUR(b);
        sub = data.pairSummary(string(data.pairSummary.UR_Band) == bn & string(data.pairSummary.Phase) == "All", :);
        sub = sortrows(sub, 'Photoperiod_h');
        col = script8_band_colour_(pal, bn);
        y = sub.Mean_Delta_log10;
        x = 1:height(sub);
        err = sub.SD_Delta_log10 ./ sqrt(max(data.nMice, 1));
        errorbar(axC, x, y, err, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, 'CapSize', 8);
        script8_direct_line_label_(axC, x(end), y(end), script8_band_display_(bn, 'tex'), col, theme, true);
        if script8_lme_pp_sig_(data.lmeFdr, bn, 'Delta')
            plot(axC, x(end), y(end), 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off');
        end
    end
    script8_yline_zero_(axC);
    set(axC, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axC, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axC, '\Delta(UR - CR) log_{10} power', 'FontWeight', 'bold');
    title(axC, 'CR–UR co-expression across photoperiod', 'FontWeight', 'bold');
    script8_style_axes_(axC, theme);
    script8_panel_label_(axC, 'C', theme);
    script8_add_n_annotation_(axC, data.nMice, [], theme);
    xl = xlim(axC); yl = ylim(axC);
    text(axC, xl(2), yl(1) + 0.08 * range(yl), 'L24: CR collapse \rightarrow UR > CR', ...
        'HorizontalAlignment', 'right', 'FontName', theme.fontName, 'FontSize', 9, 'Color', pal.l24, 'Interpreter', 'tex');
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, 'Fig03_C_Delta'), theme, {theme.ext, '.pdf'});
    close(figC);

    %% D — Sex-stratified for both primary URs (1x2 facets)
    figD = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 420]);
    for bi = 1:numel(data.primaryUR)
        axD = subplot(1, 2, bi); hold(axD, 'on');
        bn = data.primaryUR(bi);
        if ~isempty(data.pairBySex)
            pooled = data.pairSummary(string(data.pairSummary.UR_Band) == bn & string(data.pairSummary.Phase) == "All", :);
            pooled = sortrows(pooled, 'Photoperiod_h');
            plot(axD, 1:height(pooled), pooled.Mean_Delta_log10, '--', 'Color', pal.pooled, 'LineWidth', 1.2);
            for sx = ["Female", "Male"]
                sub = data.pairBySex(string(data.pairBySex.UR_Band) == bn & string(data.pairBySex.Phase) == "All" & string(data.pairBySex.Sex) == sx, :);
                sub = sortrows(sub, 'Photoperiod_h');
                if sx == "Female", col = pal.female; else, col = pal.male; end
                plot(axD, 1:height(sub), sub.Mean_Delta_log10, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col);
                script8_direct_line_label_(axD, height(sub), sub.Mean_Delta_log10(end), char(sx), col, theme, false);
            end
        end
        script8_yline_zero_(axD);
        set(axD, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
        xlabel(axD, 'Photoperiod', 'FontWeight', 'bold');
        if bi == 1
            ylabel(axD, '\Delta(UR - CR) log_{10}', 'FontWeight', 'bold');
            script8_panel_label_(axD, 'D', theme);
        end
        title(axD, ['Sex-stratified \Delta — ' script8_band_display_(bn, 'tex')], 'FontWeight', 'bold', 'Interpreter', 'tex');
        script8_style_axes_(axD, theme);
    end
    standPaths{4} = script8_export_figure_(figD, fullfile(standDir, 'Fig03_D_Sex'), theme, {theme.ext, '.pdf'});
    close(figD);

    for k = 1:3
        manifest = script8_manifest_add_(manifest, 'Fig03', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            'Gradient CR–UR panel', 'Script 6', 'Tol contract', '');
    end
    manifest = script8_manifest_add_(manifest, 'Fig03', 'D', standPaths{4}{1}, 'standalone', ...
        'Sex-stratified \Delta (Supplementary)', 'Script 6', 'Tol sex', 'Moved off main composite.');
    [widePath, tallPath] = script8_composite_fig03_(standPaths(1:3), outDirs, theme);
    manifest = script8_manifest_add_(manifest, 'Fig03', 'Composite', widePath, '16:9', ...
        'Photoperiod gradient of CR–UR co-expression with dual-band absolute power (sex in Supplementary).', ...
        'Script 6', 'Tol', 'Main A–C only.');
end

function [widePath, tallPath] = script8_composite_fig03_(standPaths, outDirs, theme)
    imgs = cell(3, 1);
    for i = 1:3
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script8_trim_image_whitespace_(imread(p));
    end
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 900]);
    tl = tiledlayout(figW, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1); imshow(imgs{1}); axis off;
    nexttile(tl, 2); imshow(imgs{2}); axis off;
    nexttile(tl, 3, [1 2]); imshow(imgs{3}); axis off;
    widePath = fullfile(outDirs.compositeWide, ['Fig03_Gradient_CRUR' theme.ext]);
    exportgraphics(figW, widePath, 'Resolution', theme.dpi);
    close(figW);

    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1400]);
    tl2 = tiledlayout(figT, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:3
        nexttile(tl2, i); imshow(imgs{i}); axis off;
    end
    tallPath = fullfile(outDirs.compositeTall, ['Fig03_Gradient_CRUR' theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
end

%% Fig04 — Transitions: C coherence hero, then A|B summary gradients
function [manifest, widePath, tallPath, standPaths] = script8_build_fig04_(data, outDirs, theme, cfg, manifest) %#ok<INUSD>
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig04_Transitions');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(3, 1);
    pp = data.ppOrder;
    ppPlot = pp(pp <= 22);

    grad = data.resyncGradient;
    if isempty(grad), grad = table(); end

    %% A — ridge power flip
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA); hold(axA, 'on');
    script8_plot_transition_metric_(axA, grad, data.primaryUR, 'RidgePower_PostMinusPre', ppPlot, pal, theme);
    xlabel(axA, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axA, 'Ridge power post - pre', 'FontWeight', 'bold');
    title(axA, 'UR ridge amplitude flip at transition', 'FontWeight', 'bold');
    script8_style_axes_(axA, theme);
    script8_panel_label_(axA, 'A', theme);
    script8_add_n_annotation_(axA, data.nMice, [], theme);
    text(axA, 0.98, 0.12, 'Stops at L22 (no entrained transition at LL)', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig04_A_RidgePowerGrad'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — phase concentration flip (\Delta R)
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    script8_plot_transition_metric_(axB, grad, data.primaryUR, 'DeltaR', ppPlot, pal, theme);
    xlabel(axB, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axB, '\DeltaR (post - pre)', 'FontWeight', 'bold');
    title(axB, 'UR phase organisation flip (\DeltaR)', 'FontWeight', 'bold');
    script8_style_axes_(axB, theme);
    script8_panel_label_(axB, 'B', theme);
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig04_B_DeltaRGrad'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — 2x3 facets L12–L22 for UR 1–3 (shared auto ylim; white faces)
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1400 780]);
    facets = pal.coherenceFacets;
    yMax = script8_coherence_ymax_(data.binnedCoherence, facets, data.primaryUR(1), pal);
    for fi = 1:numel(facets)
        axC = subplot(2, 3, fi); hold(axC, 'on');
        set(axC, 'Color', 'w');
        hasFacet = script8_plot_coherence_facet_(axC, data.binnedCoherence, facets(fi), data.primaryUR(1), pal, theme, false);
        title(axC, char(script8_pp_label_(facets(fi))), 'FontWeight', 'bold');
        if ~hasFacet
            script8_plot_empty_coherence_note_(axC, facets(fi), theme);
        end
        if fi == 1 || fi == 4
            ylabel(axC, 'Phase coherence R', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(axC, 'Time relative to transition (h)', 'FontWeight', 'bold');
        end
        if fi == 1
            script8_panel_label_(axC, 'C', theme);
        end
        set(axC, 'YLim', [0 yMax], 'XLim', pal.coherenceXlim, 'Color', 'w');
        script8_style_axes_(axC, theme);
    end
    sgtitle(figC, ['DL/LD phase coherence (' script8_band_display_(data.primaryUR(1), 'tex') ...
        '; L12–L22 entrained photoperiods)'], 'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'tex');
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, 'Fig04_C_CoherenceFacets'), theme, {theme.ext, '.pdf'});
    close(figC);

    for k = 1:3
        manifest = script8_manifest_add_(manifest, 'Fig04', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            'Transition dynamics panel', 'Script 5', 'Tol DL/LD', '');
    end
    [widePath, tallPath] = script8_composite_fig04_(standPaths, outDirs, theme);
    manifest = script8_manifest_add_(manifest, 'Fig04', 'Composite', widePath, '16:9', ...
        'Transition coherence hero (C) with ridge-power and \DeltaR summary gradients (A|B).', ...
        'Script 5', 'Tol', 'C-first layout; shared auto ylim; no grey facet fills.');
end

function [widePath, tallPath] = script8_composite_fig04_(standPaths, outDirs, theme)
    imgs = cell(3, 1);
    for i = 1:3
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script8_trim_image_whitespace_(imread(p));
    end
    % Wide: C on top full width, A|B below
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 1000]);
    tl = tiledlayout(figW, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1, [1 2]); imshow(imgs{3}); axis off;
    nexttile(tl, 3); imshow(imgs{1}); axis off;
    nexttile(tl, 4); imshow(imgs{2}); axis off;
    widePath = fullfile(outDirs.compositeWide, ['Fig04_Transitions' theme.ext]);
    exportgraphics(figW, widePath, 'Resolution', theme.dpi);
    close(figW);

    % Tall primary: C, then A, then B
    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1400]);
    tl2 = tiledlayout(figT, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl2, 1); imshow(imgs{3}); axis off;
    nexttile(tl2, 2); imshow(imgs{1}); axis off;
    nexttile(tl2, 3); imshow(imgs{2}); axis off;
    tallPath = fullfile(outDirs.compositeTall, ['Fig04_Transitions' theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
end

%% Fig05 — Cluster callouts + transition-aligned activity resync
function [manifest, widePath, tallPath, standPaths] = script8_build_fig05_(data, outDirs, theme, cfg, manifest) %#ok<INUSD>
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig05_ActivityResync');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(4, 1);

    primaryBand = data.primaryUR(1);
    primaryCluster = script8_primary_cluster_(data.clusterSummary, primaryBand);

    %% A — Cluster period callouts (lines at actual PeriodCentre_h)
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA); hold(axA, 'on');
    CS = data.clusterSummary;
    if ~isempty(CS) && all(ismember({'PeriodCentre_h', 'CandidateCount', 'BandName'}, CS.Properties.VariableNames))
        centres = double(CS.PeriodCentre_h);
        counts = double(CS.CandidateCount);
        for i = 1:numel(centres)
            bn = string(CS.BandName(i));
            col = script8_band_colour_(pal, bn);
            bar(axA, centres(i), counts(i), 0.32, 'FaceColor', col, 'EdgeColor', 'k', 'LineWidth', 1);
            xline(axA, centres(i), '--', 'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
            nMiceTxt = '';
            if ismember('MouseCount', CS.Properties.VariableNames)
                nMiceTxt = sprintf('; %d mice', double(CS.MouseCount(i)));
            end
            text(axA, centres(i), counts(i) + max(counts) * 0.04, ...
                sprintf('%s\nn=%d candidates%s', script8_band_display_(bn, 'plain'), counts(i), nMiceTxt), ...
                'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontSize', 8, 'Interpreter', 'none');
        end
        xlim(axA, [min(centres) - 0.55, max(centres) + 0.55]);
    end
    xlabel(axA, 'Validated ridge period (h)', 'FontWeight', 'bold');
    ylabel(axA, 'Validated candidates', 'FontWeight', 'bold');
    title(axA, 'Validated period clusters (UR 1–3 and UR 3–6)', 'FontWeight', 'bold');
    script8_style_axes_(axA, theme);
    script8_panel_label_(axA, 'A', theme);
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig05_A_Clusters'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — L12 transition-aligned activity (LD primary | DL secondary)
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1100 420]);
    axB1 = subplot(1, 2, 1); hold(axB1, 'on'); set(axB1, 'Color', 'w');
    script8_plot_transition_activity_(axB1, data.activityZT, primaryCluster, 12, "LD", pal, theme, true);
    title(axB1, 'L12 — LD (lights-off)', 'FontWeight', 'bold');
    ylabel(axB1, 'Activity (z-scored)', 'FontWeight', 'bold');
    script8_panel_label_(axB1, 'B', theme);
    axB2 = subplot(1, 2, 2); hold(axB2, 'on'); set(axB2, 'Color', 'w');
    script8_plot_transition_activity_(axB2, data.activityZT, primaryCluster, 12, "DL", pal, theme, true);
    title(axB2, 'L12 — DL (lights-on)', 'FontWeight', 'bold');
    xlabel(axB1, 'Time relative to transition (h)', 'FontWeight', 'bold');
    xlabel(axB2, 'Time relative to transition (h)', 'FontWeight', 'bold');
    sgtitle(figB, sprintf('Transition-aligned activity — %s (%s)', ...
        script8_band_display_(primaryBand, 'plain'), char(primaryCluster)), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig05_B_Activity_L12'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — L22 transition-aligned activity
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1100 420]);
    axC1 = subplot(1, 2, 1); hold(axC1, 'on'); set(axC1, 'Color', 'w');
    script8_plot_transition_activity_(axC1, data.activityZT, primaryCluster, 22, "LD", pal, theme, true);
    title(axC1, 'L22 — LD (lights-off)', 'FontWeight', 'bold');
    ylabel(axC1, 'Activity (z-scored)', 'FontWeight', 'bold');
    script8_panel_label_(axC1, 'C', theme);
    axC2 = subplot(1, 2, 2); hold(axC2, 'on'); set(axC2, 'Color', 'w');
    script8_plot_transition_activity_(axC2, data.activityZT, primaryCluster, 22, "DL", pal, theme, true);
    title(axC2, 'L22 — DL (lights-on)', 'FontWeight', 'bold');
    xlabel(axC1, 'Time relative to transition (h)', 'FontWeight', 'bold');
    xlabel(axC2, 'Time relative to transition (h)', 'FontWeight', 'bold');
    sgtitle(figC, sprintf('Transition-aligned activity — %s (%s)', ...
        script8_band_display_(primaryBand, 'plain'), char(primaryCluster)), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, 'Fig05_C_Activity_L22'), theme, {theme.ext, '.pdf'});
    close(figC);

    %% D — Dual-band mean overlay at L12 LD
    figD = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 720 480]);
    axD = axes(figD); hold(axD, 'on'); set(axD, 'Color', 'w');
    for bi = 1:numel(data.primaryUR)
        bn = data.primaryUR(bi);
        cid = script8_primary_cluster_(data.clusterSummary, bn);
        col = script8_band_colour_(pal, bn);
        script8_plot_transition_activity_mean_(axD, data.activityZT, cid, 12, "LD", col, theme, ...
            script8_band_display_(bn, 'tex'));
    end
    xline(axD, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.4, 'HandleVisibility', 'off');
    yline(axD, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
    xlim(axD, pal.coherenceXlim);
    xlabel(axD, 'Time relative to LD (h)', 'FontWeight', 'bold');
    ylabel(axD, 'Mean activity (z-scored)', 'FontWeight', 'bold');
    title(axD, 'Dual-band mean activity at L12 LD', 'FontWeight', 'bold');
    script8_style_axes_(axD, theme);
    script8_panel_label_(axD, 'D', theme);
    standPaths{4} = script8_export_figure_(figD, fullfile(standDir, 'Fig05_D_DualBand_L12_LD'), theme, {theme.ext, '.pdf'});
    close(figD);

    for k = 1:4
        manifest = script8_manifest_add_(manifest, 'Fig05', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            'Cluster-anchored activity resync', 'Scripts 7/8', 'Tol', '');
    end
    [widePath, tallPath] = script8_composite_from_paths_(standPaths, outDirs, 'Fig05_ActivityResync', theme, [2 2], [16 9], [4 5]);
    manifest = script8_manifest_add_(manifest, 'Fig05', 'Composite', widePath, '16:9', ...
        'Cluster-anchored transition activity resync (grey individuals, thick mean); replaces prior coherence-repeat Fig05.', ...
        'Scripts 7/8', 'Tol', 'LD primary | DL secondary; L12 vs L22.');
end

function cid = script8_primary_cluster_(CS, bandName)
    cid = "";
    if isempty(CS) || ~ismember('BandName', CS.Properties.VariableNames)
        return;
    end
    sub = CS(string(CS.BandName) == string(bandName), :);
    if isempty(sub), return; end
    if ismember('ClusterRank', sub.Properties.VariableNames)
        sub = sortrows(sub, 'ClusterRank', 'ascend');
    elseif ismember('CandidateCount', sub.Properties.VariableNames)
        sub = sortrows(sub, 'CandidateCount', 'descend');
    end
    cid = string(sub.ClusterID(1));
end

function script8_plot_transition_activity_(ax, Act, clusterID, photo, transitionType, pal, theme, showIndividuals)
    if isempty(Act) || strlength(string(clusterID)) == 0
        text(ax, 0.5, 0.5, 'No activity data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        return;
    end
    needed = {'ClusterID', 'Photoperiod_h', 'SignalID', 'ZTBinCenter_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames))
        text(ax, 0.5, 0.5, 'Activity table missing columns', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        return;
    end
    A = Act(string(Act.ClusterID) == string(clusterID) & Act.Photoperiod_h == photo, :);
    if isempty(A)
        text(ax, 0.5, 0.5, sprintf('No rows for %s / L%d', char(clusterID), photo), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'Interpreter', 'none');
        return;
    end
    transitionZT = 0;
    if upper(string(transitionType)) == "LD"
        transitionZT = double(photo);
    end
    xlimRange = theme.palette.coherenceXlim;
    meanCol = theme.palette.ld;
    if upper(string(transitionType)) == "DL"
        meanCol = theme.palette.dl;
    end
    hold(ax, 'on');
    if ~ismember('File', A.Properties.VariableNames)
        A.File = repmat("", height(A), 1);
    end
    mouseKey = string(A.File) + "|" + string(A.SignalID);
    sigs = unique(mouseKey, 'stable');
    mouseRel = {};
    mouseY = {};
    for s = 1:numel(sigs)
        As = A(mouseKey == sigs(s), :);
        As = sortrows(As, 'ZTBinCenter_h');
        zt = double(As.ZTBinCenter_h);
        y = double(As.Activity_zscored);
        rel = mod(zt - transitionZT + 12, 24) - 12;
        [rel, ord] = sort(rel);
        y = y(ord);
        keep = rel >= xlimRange(1) & rel <= xlimRange(2) & isfinite(y);
        if ~any(keep), continue; end
        if showIndividuals
            plot(ax, rel(keep), y(keep), '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 0.85, 'HandleVisibility', 'off');
        end
        mouseRel{end+1,1} = rel(keep); %#ok<AGROW>
        mouseY{end+1,1} = y(keep); %#ok<AGROW>
    end
    if ~isempty(mouseRel)
        edges = xlimRange(1):0.5:xlimRange(2);
        centers = edges(1:end-1) + diff(edges)/2;
        meanY = nan(size(centers));
        for i = 1:numel(centers)
            vals = nan(numel(mouseRel), 1);
            for m = 1:numel(mouseRel)
                idx = mouseRel{m} >= edges(i) & mouseRel{m} < edges(i+1);
                if i == numel(centers)
                    idx = mouseRel{m} >= edges(i) & mouseRel{m} <= edges(i+1);
                end
                if any(idx)
                    vals(m) = mean(mouseY{m}(idx), 'omitnan');
                end
            end
            meanY(i) = mean(vals, 'omitnan');
        end
        plot(ax, centers, meanY, '-', 'Color', meanCol, 'LineWidth', 2.6);
    end
    xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.4, 'HandleVisibility', 'off');
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
    xlim(ax, xlimRange);
    script8_style_axes_(ax, theme);
end

function script8_plot_transition_activity_mean_(ax, Act, clusterID, photo, transitionType, colour, theme, labelTxt)
    if isempty(Act) || strlength(string(clusterID)) == 0, return; end
    needed = {'ClusterID', 'Photoperiod_h', 'SignalID', 'ZTBinCenter_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames)), return; end
    A = Act(string(Act.ClusterID) == string(clusterID) & Act.Photoperiod_h == photo, :);
    if isempty(A), return; end
    transitionZT = 0;
    if upper(string(transitionType)) == "LD"
        transitionZT = double(photo);
    end
    xlimRange = theme.palette.coherenceXlim;
    if ~ismember('File', A.Properties.VariableNames)
        A.File = repmat("", height(A), 1);
    end
    mouseKey = string(A.File) + "|" + string(A.SignalID);
    sigs = unique(mouseKey, 'stable');
    mouseRel = {};
    mouseY = {};
    for s = 1:numel(sigs)
        As = A(mouseKey == sigs(s), :);
        zt = double(As.ZTBinCenter_h);
        y = double(As.Activity_zscored);
        rel = mod(zt - transitionZT + 12, 24) - 12;
        keep = isfinite(rel) & isfinite(y) & rel >= xlimRange(1) & rel <= xlimRange(2);
        if ~any(keep), continue; end
        mouseRel{end+1,1} = rel(keep); %#ok<AGROW>
        mouseY{end+1,1} = y(keep); %#ok<AGROW>
    end
    if isempty(mouseRel), return; end
    edges = xlimRange(1):0.5:xlimRange(2);
    centers = edges(1:end-1) + diff(edges)/2;
    meanY = nan(size(centers));
    for i = 1:numel(centers)
        vals = nan(numel(mouseRel), 1);
        for m = 1:numel(mouseRel)
            idx = mouseRel{m} >= edges(i) & mouseRel{m} < edges(i+1);
            if i == numel(centers)
                idx = mouseRel{m} >= edges(i) & mouseRel{m} <= edges(i+1);
            end
            if any(idx)
                vals(m) = mean(mouseY{m}(idx), 'omitnan');
            end
        end
        meanY(i) = mean(vals, 'omitnan');
    end
    plot(ax, centers, meanY, '-', 'Color', colour, 'LineWidth', 2.4);
    if any(isfinite(meanY))
        ip = find(isfinite(meanY), 1, 'last');
        script8_direct_line_label_(ax, centers(ip), meanY(ip), labelTxt, colour, theme, true);
    end
end

function yMax = script8_coherence_ymax_(Bin, facets, bandName, pal)
    yMax = pal.coherenceYMax;
    if isempty(Bin) || ~ismember('R', Bin.Properties.VariableNames)
        return;
    end
    B = Bin(ismember(Bin.Photoperiod_h, facets) & string(Bin.BandName) == string(bandName), :);
    if isempty(B), return; end
    mx = max(double(B.R), [], 'omitnan');
    if isfinite(mx) && mx > 0
        yMax = max(ceil(mx * 1.05 * 20) / 20, 0.2); % neat ceiling ~0.05 steps
    end
end

function script8_plot_binned_ridge_(ax, Bin, bandName, photos, pal, theme)
    if isempty(Bin), return; end
    hold(ax, 'on');
    cols = {pal.l12, pal.l24};
    for i = 1:numel(photos)
        B = Bin(Bin.Photoperiod_h == photos(i) & string(Bin.BandName) == bandName & ...
            (string(Bin.TransitionType) == "LD" | string(Bin.TransitionType) == "DL"), :);
        if isempty(B) || ~ismember('MeanRidgePower_log10', B.Properties.VariableNames)
            continue;
        end
        % Average LD and DL within relative-time bins for a single photoperiod trace
        B = sortrows(B, 'RelBinCenter_h');
        G = groupsummary(B, 'RelBinCenter_h', 'mean', 'MeanRidgePower_log10');
        mcol = "mean_MeanRidgePower_log10";
        if ~ismember(mcol, G.Properties.VariableNames)
            vn = G.Properties.VariableNames(startsWith(G.Properties.VariableNames, 'mean_'));
            if isempty(vn), continue; end
            mcol = string(vn{1});
        end
        plot(ax, G.RelBinCenter_h, G.(mcol), '-', 'Color', cols{i}, 'LineWidth', 2.0);
        script8_direct_line_label_(ax, G.RelBinCenter_h(end), G.(mcol)(end), ...
            char(script8_pp_label_(photos(i))), cols{i}, theme, false);
    end
    xline(ax, 0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlim(ax, theme.palette.coherenceXlim);
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
            script8_direct_line_label_(ax, ip, vals(ip), char(tr), col, theme, false);
        end
    end
end

function hasData = script8_plot_coherence_facet_(ax, Bin, photoperiod_h, bandName, pal, theme, doShade)
    if nargin < 7, doShade = false; end
    hasData = false;
    if isempty(Bin), return; end
    B = Bin(Bin.Photoperiod_h == photoperiod_h & string(Bin.BandName) == bandName, :);
    if isempty(B), return; end
    hold(ax, 'on');
    set(ax, 'Color', 'w');
    for tr = ["DL", "LD"]
        Bt = B(string(B.TransitionType) == tr, :);
        if isempty(Bt), continue; end
        hasData = true;
        Bt = sortrows(Bt, 'RelBinCenter_h');
        col = script8_transition_colour_(pal, tr);
        x = Bt.RelBinCenter_h;
        y = Bt.R;
        plot(ax, x, y, '-o', 'Color', col, 'LineWidth', 2.0, 'MarkerSize', 4, 'MarkerFaceColor', col);
        if ismember('N_PhaseObs', Bt.Properties.VariableNames)
            se = 1.96 ./ sqrt(max(Bt.N_PhaseObs, 1));
            fill(ax, [x; flipud(x)], [y - se; flipud(y + se)], col, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        end
        if x(end) > 3
            script8_direct_line_label_(ax, x(end), y(end), char(tr), col, theme, false);
        end
        if doShade
            script8_shade_transition_(ax, pal.coherenceXlim, tr == "DL", pal);
        end
    end
    xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.3, 'HandleVisibility', 'off');
    script8_style_axes_(ax, theme);
end

function script8_plot_empty_coherence_note_(ax, photoperiod_h, theme)
    axis(ax, [0 1 0 1]); axis(ax, 'off');
    text(ax, 0.5, 0.55, sprintf('No data for %s', char(script8_pp_label_(photoperiod_h))), ...
        'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 10);
end

function script8_plot_24h_profile_(ax, T, yCol, photos, pal, theme, yLabel, bandName)
    if nargin < 8, bandName = "UR_1_3"; end
    if isempty(T) || ~ismember(yCol, T.Properties.VariableNames)
        return;
    end
    if ismember('BandName', T.Properties.VariableNames)
        T = T(string(T.BandName) == bandName, :);
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
        labelX = subG.(xCol)(max(numel(subG.(xCol)) - 1, 1));
        labelY = subG.(meanCol)(max(numel(subG.(meanCol)) - 1, 1));
        script8_direct_line_label_(ax, labelX + 0.15, labelY, char(labs{i}), cols{i}, theme, false);
    end
    xlabel(ax, 'ZT (h)', 'FontWeight', 'bold');
    ylabel(ax, yLabel, 'FontWeight', 'bold');
    xlim(ax, [0 24.75]);
    script8_style_axes_(ax, theme);
end

%% Supplementary — sex panel + 24h profiles + coexpression (+ optional legacy)
function [manifest, standPaths] = script8_build_supplementary_(data, paths, outDirs, theme, manifest)
    standDir = fullfile(outDirs.standalone, 'Supplementary');
    extended_period_gate_ensure_dir(standDir);
    standPaths = {};
    pal = theme.palette;
    pp = data.ppOrder;
    ppLabels = arrayfun(@(x) char(script8_pp_label_(x)), pp, 'UniformOutput', false);

    % Sex-stratified delta (former Fig03D)
    figSex = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 420]);
    for bi = 1:numel(data.primaryUR)
        axD = subplot(1, 2, bi); hold(axD, 'on');
        bn = data.primaryUR(bi);
        if ~isempty(data.pairBySex)
            pooled = data.pairSummary(string(data.pairSummary.UR_Band) == bn & string(data.pairSummary.Phase) == "All", :);
            pooled = sortrows(pooled, 'Photoperiod_h');
            plot(axD, 1:height(pooled), pooled.Mean_Delta_log10, '--', 'Color', pal.pooled, 'LineWidth', 1.2);
            for sx = ["Female", "Male"]
                sub = data.pairBySex(string(data.pairBySex.UR_Band) == bn & string(data.pairBySex.Phase) == "All" & string(data.pairBySex.Sex) == sx, :);
                sub = sortrows(sub, 'Photoperiod_h');
                if sx == "Female", col = pal.female; else, col = pal.male; end
                plot(axD, 1:height(sub), sub.Mean_Delta_log10, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col);
                script8_direct_line_label_(axD, height(sub), sub.Mean_Delta_log10(end), char(sx), col, theme, false);
            end
        end
        script8_yline_zero_(axD);
        set(axD, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
        xlabel(axD, 'Photoperiod', 'FontWeight', 'bold');
        if bi == 1
            ylabel(axD, '\Delta(UR - CR) log_{10}', 'FontWeight', 'bold');
        end
        title(axD, ['Sex-stratified \Delta — ' script8_band_display_(bn, 'tex')], 'FontWeight', 'bold', 'Interpreter', 'tex');
        script8_style_axes_(axD, theme);
    end
    outSex = script8_export_figure_(figSex, fullfile(standDir, 'Supp_Fig03D_Sex'), theme, {theme.ext, '.pdf'});
    close(figSex);
    standPaths{end + 1} = outSex{1}; %#ok<AGROW>
    manifest = script8_manifest_add_(manifest, 'Supp', 'Sex_delta', outSex{1}, 'standalone', ...
        'Sex-stratified CR–UR \Delta (moved from main Fig03).', 'Script 6', 'Tol sex', '');

    % 24h L12 vs L24 profiles
    figP = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 420]);
    ax1 = subplot(1, 2, 1); hold(ax1, 'on');
    script8_plot_24h_profile_(ax1, data.phase24, 'R', [12 24], pal, theme, 'Phase coherence R', "UR_1_3");
    title(ax1, ['24h phase coherence — ' script8_band_display_('UR_1_3', 'tex')], 'FontWeight', 'bold', 'Interpreter', 'tex');
    ax2 = subplot(1, 2, 2); hold(ax2, 'on');
    script8_plot_24h_profile_(ax2, data.ridge24, 'MeanRidgePower_log10', [12 24], pal, theme, 'Ridge power (log_{10})', "UR_1_3");
    title(ax2, ['24h ridge power — ' script8_band_display_('UR_1_3', 'tex')], 'FontWeight', 'bold', 'Interpreter', 'tex');
    outP = script8_export_figure_(figP, fullfile(standDir, 'Supp_24h_L12_vs_L24'), theme, {theme.ext, '.pdf'});
    close(figP);
    standPaths{end + 1} = outP{1}; %#ok<AGROW>
    manifest = script8_manifest_add_(manifest, 'Supp', '24h_profiles', outP{1}, 'standalone', ...
        '24h L12 vs L24 phase/ridge profiles (UR 1–3).', 'Script 7', 'Tol', 'Moved off main Fig05.');

    % Legacy dual-band coherence / ridge-at-transition (former Fig05 B–D)
    figLegacy = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1200 780]);
    heroPP = [12 16 22];
    for fi = 1:numel(heroPP)
        axB = subplot(2, 3, fi); hold(axB, 'on'); set(axB, 'Color', 'w');
        script8_plot_coherence_facet_(axB, data.binnedCoherence, heroPP(fi), "UR_1_3", pal, theme, false);
        title(axB, [char(script8_pp_label_(heroPP(fi))) ' UR 1–3'], 'FontWeight', 'bold');
        set(axB, 'YLim', [0 pal.coherenceYMax], 'XLim', pal.coherenceXlim);
        script8_style_axes_(axB, theme);
        axC = subplot(2, 3, fi + 3); hold(axC, 'on'); set(axC, 'Color', 'w');
        script8_plot_coherence_facet_(axC, data.binnedCoherence, heroPP(fi), "UR_3_6", pal, theme, false);
        title(axC, [char(script8_pp_label_(heroPP(fi))) ' UR 3–6'], 'FontWeight', 'bold');
        set(axC, 'YLim', [0 pal.coherenceYMax], 'XLim', pal.coherenceXlim);
        script8_style_axes_(axC, theme);
    end
    outL = script8_export_figure_(figLegacy, fullfile(standDir, 'Supp_Legacy_DualBand_Coherence'), theme, {theme.ext, '.pdf'});
    close(figLegacy);
    standPaths{end + 1} = outL{1}; %#ok<AGROW>
    manifest = script8_manifest_add_(manifest, 'Supp', 'Legacy_coherence', outL{1}, 'standalone', ...
        'Former Fig05 dual-band coherence facets (now Supplementary).', 'Scripts 5/7', 'Tol', '');

    figRidge = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 420]);
    for bi = 1:2
        axD = subplot(1, 2, bi); hold(axD, 'on');
        bn = data.primaryUR(bi);
        script8_plot_binned_ridge_(axD, data.binnedCoherence, bn, [12 22], pal, theme);
        title(axD, ['Ridge power at LD/DL — ' script8_band_display_(bn, 'tex')], ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
        if bi == 1
            ylabel(axD, 'Mean ridge power (log_{10})', 'FontWeight', 'bold');
        end
        xlabel(axD, 'Time relative to transition (h)', 'FontWeight', 'bold');
        script8_style_axes_(axD, theme);
    end
    outR = script8_export_figure_(figRidge, fullfile(standDir, 'Supp_Legacy_RidgeAtTransition'), theme, {theme.ext, '.pdf'});
    close(figRidge);
    standPaths{end + 1} = outR{1}; %#ok<AGROW>
    manifest = script8_manifest_add_(manifest, 'Supp', 'Legacy_ridge', outR{1}, 'standalone', ...
        'Former Fig05 ridge-at-transition panels (now Supplementary).', 'Scripts 5/7', 'Tol', '');

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
        if isfile(p), imgs{i} = script8_trim_image_whitespace_(imread(p)); end
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
        if iscell(img), img = script8_trim_image_whitespace_(imread(img{1})); end
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
        if iscell(img), img = script8_trim_image_whitespace_(imread(img{1})); end
        if ~isempty(img), imshow(img); axis off; end
    end
    tallPath = fullfile(outDirs.compositeTall, [name theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
end

function script8_export_tall_copy_(tallPath, widePath, tallDir)
    if isfile(tallPath), return; end
    if isfile(widePath)
        copyfile(widePath, fullfile(tallDir, ['TALL_' basename_(widePath)]), 'f');
    end
end

function b = basename_(p)
    [~, b, e] = fileparts(p);
    b = [b e];
end

function img = script8_trim_image_whitespace_(img)
    if isempty(img), return; end
    if ndims(img) == 3
        mask = any(img < 250, 3);
    else
        mask = img < 250;
    end
    [rows, cols] = find(mask);
    if isempty(rows) || isempty(cols), return; end
    pad = 6;
    r1 = max(min(rows) - pad, 1);
    r2 = min(max(rows) + pad, size(img, 1));
    c1 = max(min(cols) - pad, 1);
    c2 = min(max(cols) + pad, size(img, 2));
    img = img(r1:r2, c1:c2, :);
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
        'UR 1–3 h', '#66CCEE', 'Band panels (key UR_1_3)';
        'UR 3–6 h', '#AA3377', 'Band panels (key UR_3_6)';
        'CR 20–28 h', '#BBBBBB', 'CR traces / dark schematic';
        'Female', '#228833', 'Sex inset only';
        'Male', '#CCBB44', 'Sex inset only';
        'L12 (compare)', '#4477AA', 'Fig05 / supp 24h';
        'L22 / L24', '#AA3377', 'Fig05 / supp 24h';
        'Scalograms', 'Jet', 'Fig02 RAW + HSub residual'};
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

function lbl = script8_band_display_(bandKey, mode)
%SCRIPT8_BAND_DISPLAY_ Map machine keys to figure-facing labels.
    if nargin < 2, mode = 'tex'; end
    key = char(string(bandKey));
    switch key
        case 'UR_1_3'
            texLbl = 'UR_{1–3}'; plainLbl = 'UR 1–3 h';
        case 'UR_3_6'
            texLbl = 'UR_{3–6}'; plainLbl = 'UR 3–6 h';
        case 'UR_6_9'
            texLbl = 'UR_{6–9}'; plainLbl = 'UR 6–9 h';
        case 'UR_9_12'
            texLbl = 'UR_{9–12}'; plainLbl = 'UR 9–12 h';
        case 'UR_12_18'
            texLbl = 'UR_{12–18}'; plainLbl = 'UR 12–18 h';
        case 'CR_20_28'
            texLbl = 'CR_{20–28}'; plainLbl = 'CR 20–28 h';
        otherwise
            texLbl = strrep(key, '_', '\_'); plainLbl = strrep(key, '_', ' ');
    end
    if strcmpi(mode, 'plain')
        lbl = plainLbl;
    else
        lbl = texLbl;
    end
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

function script8_direct_line_label_(ax, x, y, txt, colour, theme, useTex)
    if nargin < 7, useTex = false; end
    interp = 'none';
    if useTex, interp = 'tex'; end
    text(ax, x, y, txt, 'Color', colour, 'FontName', theme.fontName, 'FontWeight', 'bold', ...
        'FontSize', 10, 'Interpreter', interp);
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

function script8_axes_arrow_(ax, x1, y1, x2, y2, col)
% Draw arrow in axes data coordinates (avoids rogue figure-annotation arrows).
    hold(ax, 'on');
    dx = x2 - x1; dy = y2 - y1;
    segLen = hypot(dx, dy);
    if segLen < 1e-6, return; end
    ang = atan2(dy, dx);
    headLen = min(0.10, 0.40 * segLen);
    headW = min(0.045, 0.55 * headLen);
    xBase = x2 - headLen * cos(ang);
    yBase = y2 - headLen * sin(ang);
    plot(ax, [x1 xBase], [y1 yBase], '-', 'Color', col, 'LineWidth', 1.6, 'HandleVisibility', 'off');
    px = [-headW * sin(ang), headW * sin(ang)];
    py = [headW * cos(ang), -headW * cos(ang)];
    patch(ax, [x2, xBase + px(1), xBase + px(2)], [y2, yBase + py(1), yBase + py(2)], col, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
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
