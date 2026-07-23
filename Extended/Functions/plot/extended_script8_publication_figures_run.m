function out = extended_script8_publication_figures_run(cohortRoot, cfg)
%EXTENDED_SCRIPT8_PUBLICATION_FIGURES_RUN Curated publication composites (Scripts 1-7 inputs).
%
%   Main figures:
%     Fig01a High-level methods schematic (5 steps)
%     Fig01b Detailed methods / validation schematic
%     Fig02 RAW + HSub scalograms + CarryForward retention (A-E)
%     Fig03 Photoperiod gradient CR–UR (A-C; sex in Supplementary)
%     Fig04 UR 1–3 transitions (A coherence, B amplitude, C DeltaR)
%     Fig05 UR 3–6 transitions (same layout as Fig04)
%     Fig06 Clusters + UR 1–3 primary-cluster 24h activity L12–L22
%   Supplementary: cluster identity cards (UR 1–3 + UR 3–6 primary clusters),
%     sex-stratified delta, non-primary activity grids, legacy coherence/ridge, coexpression

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

    %% Fig01a — High-level methods
    [manifest, wide1a, ~] = script8_build_fig01a_(paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide1a); %#ok<AGROW>
    script8_log_(LOG, 'Fig 01a complete');

    %% Fig01b — Detailed methods
    [manifest, wide1b, ~] = script8_build_fig01b_(paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide1b); %#ok<AGROW>
    script8_log_(LOG, 'Fig 01b complete');

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

    %% Fig04 — UR 1–3 transitions (A coherence, B|C summaries)
    band13 = data.primaryUR(1);
    [manifest, wide4, tall4] = script8_build_transition_figure_(data, outDirs, theme, cfg, manifest, ...
        'Fig04', band13, 'Fig04_Transitions_UR13');
    compositeWide(end + 1, 1) = string(wide4); %#ok<AGROW>
    script8_export_tall_copy_(tall4, wide4, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 04 complete');

    %% Fig05 — UR 3–6 transitions (twin of Fig04)
    band36 = data.primaryUR(min(2, numel(data.primaryUR)));
    [manifest, wide5, tall5] = script8_build_transition_figure_(data, outDirs, theme, cfg, manifest, ...
        'Fig05', band36, 'Fig05_Transitions_UR36');
    compositeWide(end + 1, 1) = string(wide5); %#ok<AGROW>
    script8_export_tall_copy_(tall5, wide5, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 05 complete');

    %% Fig06 — Clusters + UR 1–3 primary-cluster activity L12–L22
    [manifest, wide6, tall6] = script8_build_fig06_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide6); %#ok<AGROW>
    script8_export_tall_copy_(tall6, wide6, outDirs.compositeTall);
    script8_log_(LOG, 'Fig 06 complete');

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

%% Fig01a — High-level methods schematic (publication narrative)
function [manifest, widePath, standPaths] = script8_build_fig01a_(paths, outDirs, theme, manifest) %#ok<INUSD>
    pal = theme.palette;
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 1400 420], 'Visible', 'off');
    ax = axes(fig, 'Position', [0.04 0.22 0.92 0.66]); hold(ax, 'on');
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

    base = fullfile(outDirs.standalone, 'Fig01a_MethodsOverview');
    standPaths = script8_export_figure_(fig, base, theme, {theme.ext, '.pdf'});
    widePath = fullfile(outDirs.compositeWide, ['Fig01a_MethodsOverview' theme.ext]);
    copyfile(standPaths{1}, widePath, 'f');
    close(fig);
    manifest = script8_manifest_add_(manifest, 'Fig01a', 'All', widePath, '16:9', ...
        'High-level analysis steps: raw activity, wavelet, harmonic subtraction, CR–UR co-expression, ultradian phase coherence.', ...
        'Script8', 'Schematic', 'Script numbers omitted from boxes; see caption.');
end

%% Fig01b — Detailed methods / validation schematic
function [manifest, widePath, standPaths] = script8_build_fig01b_(paths, outDirs, theme, manifest) %#ok<INUSD>
    pal = theme.palette;
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 1500 560], 'Visible', 'off');
    ax = axes(fig, 'Position', [0.03 0.20 0.94 0.70]); hold(ax, 'on');
    axis(ax, [0 12 0 1]); axis(ax, 'off');

    w = 1.65; h = 0.30;
    yMain = 0.58;
    yBranch = 0.28;

    boxes = {
        0.9,  yMain, 'Wavelet\nscalogram';
        3.0,  yMain, 'HSub residual\n(SEL\_P360)';
        5.4,  yMain, 'CarryForward\ngate';
        7.8,  yMain, 'CR–UR\nco-expression';
        10.2, yMain, 'Transition\nphase coherence';
        };
    for i = 1:size(boxes, 1)
        x = boxes{i, 1}; y = boxes{i, 2}; txt = boxes{i, 3};
        rectangle(ax, 'Position', [x - w/2, y - h/2, w, h], 'Curvature', 0.08, ...
            'FaceColor', [1 1 1], 'EdgeColor', pal.base(1, :), 'LineWidth', 1.6);
        text(ax, x, y, txt, 'HorizontalAlignment', 'center', 'FontName', theme.fontName, ...
            'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'tex');
        if i < size(boxes, 1)
            script8_axes_arrow_(ax, x + w/2 + 0.05, y, boxes{i+1, 1} - w/2 - 0.05, boxes{i+1, 2}, pal.base(1, :));
        end
    end

    % Branch: transition metrics + clustering
    branchBoxes = {
        10.2, yBranch, '\DeltaR + ridge\nat transition';
        7.8,  yBranch, 'Validated UR\nclusters';
        5.4,  yBranch, '24h activity\ncomponents';
        };
    for i = 1:size(branchBoxes, 1)
        x = branchBoxes{i, 1}; y = branchBoxes{i, 2}; txt = branchBoxes{i, 3};
        rectangle(ax, 'Position', [x - w/2, y - h/2, w, h], 'Curvature', 0.08, ...
            'FaceColor', [1 1 1], 'EdgeColor', pal.base(3, :), 'LineWidth', 1.4);
        text(ax, x, y, txt, 'HorizontalAlignment', 'center', 'FontName', theme.fontName, ...
            'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'tex');
        if i < size(branchBoxes, 1)
            script8_axes_arrow_(ax, x - w/2 - 0.05, y, branchBoxes{i+1, 1} + w/2 + 0.05, branchBoxes{i+1, 2}, pal.base(3, :));
        end
    end
    script8_axes_arrow_(ax, 10.2, yMain - h/2 - 0.02, 10.2, yBranch + h/2 + 0.02, pal.base(3, :));

    % Wavelet → HSub link (Script 2 → Script 1)
    text(ax, 1.95, yMain + 0.14, 'Scripts 2 \rightarrow 1', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35], 'Interpreter', 'tex');
    text(ax, 4.2, yMain + 0.14, 'Script 4', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    text(ax, 6.6, yMain + 0.14, 'Scripts 3/6', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    text(ax, 9.0, yMain + 0.14, 'Script 5', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    text(ax, 6.6, yBranch - 0.18, 'Script 7', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);

    text(ax, 6, 0.06, ...
        ['CarryForward: Raw UR ridge candidates retained only when period-matched to HSub residual ' ...
        '(\pm15% SEL\_P360).'], ...
        'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontSize', 10, ...
        'FontWeight', 'bold', 'Interpreter', 'tex');
    title(ax, 'Detailed validation pipeline', 'FontName', theme.fontName, 'FontWeight', 'bold');

    base = fullfile(outDirs.standalone, 'Fig01b_MethodsDetailed');
    standPaths = script8_export_figure_(fig, base, theme, {theme.ext, '.pdf'});
    widePath = fullfile(outDirs.compositeWide, ['Fig01b_MethodsDetailed' theme.ext]);
    copyfile(standPaths{1}, widePath, 'f');
    close(fig);
    manifest = script8_manifest_add_(manifest, 'Fig01b', 'All', widePath, '16:9', ...
        'Detailed validation pipeline: wavelet scalogram, HSub residual, CarryForward gate, CR–UR co-expression, transition phase coherence, validated clusters, and 24h activity components.', ...
        'Scripts 1–7', 'Schematic', 'Scientific operation labels; script refs as footnotes.');
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

%% Fig04/Fig05 — Transitions: A coherence hero, B|C summary gradients
function [manifest, widePath, tallPath, standPaths] = script8_build_transition_figure_(data, outDirs, theme, cfg, manifest, figId, bandName, compositeStem) %#ok<INUSD>
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, compositeStem);
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(3, 1);
    pp = data.ppOrder;
    ppPlot = pp(pp <= 22);
    bandName = string(bandName);

    grad = data.resyncGradient;
    if isempty(grad), grad = table(); end

    %% A — 2x3 facets L12–L22 (shared auto ylim; white faces)
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1600 880]);
    facets = pal.coherenceFacets;
    yMax = script8_coherence_ymax_(data.binnedCoherence, facets, bandName, pal);
    for fi = 1:numel(facets)
        axA = subplot(2, 3, fi); hold(axA, 'on');
        set(axA, 'Color', 'w');
        hasFacet = script8_plot_coherence_facet_(axA, data.binnedCoherence, facets(fi), bandName, pal, theme, false);
        title(axA, char(script8_pp_label_(facets(fi))), 'FontWeight', 'bold');
        if ~hasFacet
            script8_plot_empty_coherence_note_(axA, facets(fi), theme);
        end
        if fi == 1 || fi == 4
            ylabel(axA, 'Phase coherence R', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(axA, 'Time relative to transition (h)', 'FontWeight', 'bold');
        end
        if fi == 1
            script8_panel_label_(axA, 'A', theme);
        end
        set(axA, 'YLim', [0 yMax], 'XLim', pal.coherenceXlim, 'Color', 'w');
        script8_style_axes_(axA, theme);
    end
    sgtitle(figA, ['DL/LD phase coherence (' script8_band_display_(bandName, 'tex') ...
        '; L12–L22 entrained photoperiods)'], 'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'tex');
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, [compositeStem '_A_CoherenceFacets']), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — ridge power flip
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    script8_plot_transition_metric_(axB, grad, bandName, 'RidgePower_PostMinusPre', ppPlot, pal, theme);
    xlabel(axB, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axB, 'Ridge power post - pre', 'FontWeight', 'bold');
    title(axB, ['UR ridge amplitude flip at transition — ' script8_band_display_(bandName, 'tex')], ...
        'FontWeight', 'bold', 'Interpreter', 'tex');
    script8_style_axes_(axB, theme);
    script8_panel_label_(axB, 'B', theme);
    script8_add_n_annotation_(axB, data.nMice, [], theme);
    text(axB, 0.98, 0.12, 'Stops at L22 (no entrained transition at LL)', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, [compositeStem '_B_RidgePowerGrad']), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — phase concentration flip (\Delta R)
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axC = axes(figC); hold(axC, 'on');
    script8_plot_transition_metric_(axC, grad, bandName, 'DeltaR', ppPlot, pal, theme);
    xlabel(axC, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axC, '\DeltaR (post - pre)', 'FontWeight', 'bold');
    title(axC, ['UR phase organisation flip (\DeltaR) — ' script8_band_display_(bandName, 'tex')], ...
        'FontWeight', 'bold', 'Interpreter', 'tex');
    script8_style_axes_(axC, theme);
    script8_panel_label_(axC, 'C', theme);
    standPaths{3} = script8_export_figure_(figC, fullfile(standDir, [compositeStem '_C_DeltaRGrad']), theme, {theme.ext, '.pdf'});
    close(figC);

    bandCaption = script8_band_display_(bandName, 'plain');
    for k = 1:3
        manifest = script8_manifest_add_(manifest, figId, char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            sprintf('Transition dynamics panel (%s)', bandCaption), 'Script 5', 'Tol DL/LD', '');
    end
    [widePath, tallPath] = script8_composite_transition_figure_(standPaths, outDirs, compositeStem, theme);
    manifest = script8_manifest_add_(manifest, figId, 'Composite', widePath, '16:9', ...
        sprintf('Transition coherence hero (A) with ridge-power and \\DeltaR summary gradients (B|C) for %s.', bandCaption), ...
        'Script 5', 'Tol', 'A-first layout; shared auto ylim; no grey facet fills.');
end

function [widePath, tallPath] = script8_composite_transition_figure_(standPaths, outDirs, compositeStem, theme)
    imgs = cell(3, 1);
    for i = 1:3
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script8_trim_image_whitespace_(imread(p));
    end
    % Wide: A on top full width (2 rows), B|C below
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 1100]);
    tl = tiledlayout(figW, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1, [2 2]); imshow(imgs{1}); axis off;
    nexttile(tl, 5); imshow(imgs{2}); axis off;
    nexttile(tl, 6); imshow(imgs{3}); axis off;
    widePath = fullfile(outDirs.compositeWide, [compositeStem theme.ext]);
    exportgraphics(figW, widePath, 'Resolution', theme.dpi);
    close(figW);

    % Tall: A larger, then B|C below
    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1600]);
    tl2 = tiledlayout(figT, 4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl2, 1, [3 2]); imshow(imgs{1}); axis off;
    nexttile(tl2, 7); imshow(imgs{2}); axis off;
    nexttile(tl2, 8); imshow(imgs{3}); axis off;
    tallPath = fullfile(outDirs.compositeTall, [compositeStem theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
end

%% Fig06 — Cluster callouts + primary UR 1–3 24h activity L12–L22
function [manifest, widePath, tallPath, standPaths] = script8_build_fig06_(data, outDirs, theme, cfg, manifest) %#ok<INUSD>
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig06_Clusters_Activity');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(2, 1);

    primaryBand = data.primaryUR(1);
    primaryCluster = script8_primary_cluster_(data.clusterSummary, primaryBand);
    facets = pal.coherenceFacets;
    yMaxAct = script8_activity_zt_ymax_(data.activityZT, primaryCluster, facets);

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
    standPaths{1} = script8_export_figure_(figA, fullfile(standDir, 'Fig06_A_Clusters'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — ZT 0–24 activity grid (primary UR 1–3 cluster)
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1600 880]);
    for fi = 1:numel(facets)
        axB = subplot(2, 3, fi); hold(axB, 'on');
        set(axB, 'Color', 'w');
        hasData = script8_plot_zt_activity_facet_(axB, data.activityZT, primaryCluster, facets(fi), pal, theme, yMaxAct);
        title(axB, char(script8_pp_label_(facets(fi))), 'FontWeight', 'bold');
        if ~hasData
            text(axB, 0.5, 0.5, sprintf('No activity for %s', char(script8_pp_label_(facets(fi)))), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        end
        if fi == 1 || fi == 4
            ylabel(axB, 'Activity (z-scored)', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(axB, 'ZT (h)', 'FontWeight', 'bold');
        end
        if fi == 1
            script8_panel_label_(axB, 'B', theme);
        end
        set(axB, 'XLim', [0 24], 'Color', 'w');
        script8_style_axes_(axB, theme);
    end
    sgtitle(figB, sprintf('24h activity — %s (%s)', ...
        script8_band_display_(primaryBand, 'plain'), char(primaryCluster)), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
    standPaths{2} = script8_export_figure_(figB, fullfile(standDir, 'Fig06_B_Activity_L12_L22'), theme, {theme.ext, '.pdf'});
    close(figB);

    manifest = script8_manifest_add_(manifest, 'Fig06', 'A', standPaths{1}{1}, 'standalone', ...
        'Validated period cluster callouts', 'Script 7', 'Tol bands', '');
    manifest = script8_manifest_add_(manifest, 'Fig06', 'B', standPaths{2}{1}, 'standalone', ...
        sprintf('24h ZT activity L12–L22 for primary %s cluster (%s)', ...
        script8_band_display_(primaryBand, 'plain'), script8_cluster_short_label_(primaryCluster, data.clusterSummary)), ...
        'Script 7', 'Tol L12/L22', 'Grey individuals + thick mean; LD shading; shared ylim across facets.');
    [widePath, tallPath] = script8_composite_fig06_(standPaths, outDirs, theme);
    manifest = script8_manifest_add_(manifest, 'Fig06', 'Composite', widePath, '16:9', ...
        sprintf('Validated clusters (A) and 24h activity L12–L22 for primary %s cluster (%s, panel B).', ...
        script8_band_display_(primaryBand, 'plain'), script8_cluster_short_label_(primaryCluster, data.clusterSummary)), ...
        'Script 7', 'Tol', 'Activity grid is hero panel; cluster callout compact above.');
end

function [widePath, tallPath] = script8_composite_fig06_(standPaths, outDirs, theme)
    imgs = cell(2, 1);
    for i = 1:2
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script8_trim_image_whitespace_(imread(p));
    end
    stem = 'Fig06_Clusters_Activity_UR13';
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 1280]);
    tl = tiledlayout(figW, 6, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1, [1 2]); imshow(imgs{1}); axis off;
    nexttile(tl, 3, [4 2]); imshow(imgs{2}); axis off;
    widePath = fullfile(outDirs.compositeWide, [stem theme.ext]);
    exportgraphics(figW, widePath, 'Resolution', theme.dpi);
    close(figW);

    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1680]);
    tl2 = tiledlayout(figT, 6, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl2, 1); imshow(imgs{1}); axis off;
    nexttile(tl2, 2, [5 1]); imshow(imgs{2}); axis off;
    tallPath = fullfile(outDirs.compositeTall, [stem theme.ext]);
    exportgraphics(figT, tallPath, 'Resolution', theme.dpi);
    close(figT);
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

function lbl = script8_cluster_short_label_(clusterID, CS)
    lbl = char(string(clusterID));
    if isempty(CS) || strlength(string(clusterID)) == 0
        return;
    end
    if ~ismember('ClusterID', CS.Properties.VariableNames)
        return;
    end
    sub = CS(string(CS.ClusterID) == string(clusterID), :);
    if isempty(sub)
        return;
    end
    if ismember('ClusterRank', sub.Properties.VariableNames)
        lbl = sprintf('C%02d', double(sub.ClusterRank(1)));
    end
end

function yMax = script8_coherence_zt_ymax_(Phase24, clusterIDs, photos, pal)
    yMax = [0 1];
    if nargin < 4 || isempty(pal)
        yMax = [0 0.65];
    end
    if isempty(Phase24) || ~ismember('R', Phase24.Properties.VariableNames)
        if nargin >= 4 && ~isempty(pal)
            yMax = [0 pal.coherenceYMax];
        end
        return;
    end
    clusterIDs = string(clusterIDs);
    clusterIDs = clusterIDs(strlength(clusterIDs) > 0);
    if isempty(clusterIDs)
        if nargin >= 4 && ~isempty(pal)
            yMax = [0 pal.coherenceYMax];
        end
        return;
    end
    P = Phase24(ismember(Phase24.Photoperiod_h, photos), :);
    if ismember('ClusterID', P.Properties.VariableNames)
        P = P(ismember(string(P.ClusterID), clusterIDs), :);
    end
    if isempty(P)
        if nargin >= 4 && ~isempty(pal)
            yMax = [0 pal.coherenceYMax];
        end
        return;
    end
    mx = max(double(P.R), [], 'omitnan');
    if ~isfinite(mx) || mx <= 0
        ceiling = 0.65;
        if nargin >= 4 && ~isempty(pal), ceiling = pal.coherenceYMax; end
        yMax = [0 ceiling];
        return;
    end
    yTop = max(ceil(mx * 1.05 * 20) / 20, 0.2);
    yMax = [0 yTop];
end

function [hasData, filterNote] = script8_plot_zt_coherence_facet_(ax, Phase24, clusterID, photo, lineCol, pal, theme, yLim, bandName)
    if nargin < 8, yLim = [0 1]; end
    if nargin < 9, bandName = ""; end
    hasData = false;
    filterNote = "";
    if isempty(Phase24) || strlength(string(clusterID)) == 0
        return;
    end
    needed = {'Photoperiod_h', 'R'};
    if ~all(ismember(needed, Phase24.Properties.VariableNames))
        return;
    end
    xCol = 'ZTBinCenter_h';
    if ~ismember(xCol, Phase24.Properties.VariableNames)
        if ismember('RelBinCenter_h', Phase24.Properties.VariableNames)
            xCol = 'RelBinCenter_h';
        else
            return;
        end
    end
    P = Phase24(Phase24.Photoperiod_h == photo, :);
    if ismember('ClusterID', P.Properties.VariableNames)
        P = P(string(P.ClusterID) == string(clusterID), :);
    elseif ismember('BandName', P.Properties.VariableNames) && strlength(string(bandName)) > 0
        P = P(string(P.BandName) == string(bandName), :);
        filterNote = "band-pooled (no ClusterID column)";
    end
    if isempty(P)
        return;
    end
    P = sortrows(P, xCol);
    x = double(P.(xCol));
    y = double(P.R);
    keep = x >= 0 & x <= 24 & isfinite(y);
    if ~any(keep)
        return;
    end
    hold(ax, 'on');
    set(ax, 'Color', 'w');
    plot(ax, x(keep), y(keep), '-o', 'Color', lineCol, 'LineWidth', 2.2, ...
        'MarkerSize', 4, 'MarkerFaceColor', lineCol);
    hasData = true;
    if ~isempty(yLim)
        ylim(ax, yLim);
        script8_shade_zt_ld_(ax, double(photo), pal, yLim);
    else
        yl = ylim(ax);
        script8_shade_zt_ld_(ax, double(photo), pal, yl);
    end
    photoH = double(photo);
    xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    if photoH > 0 && photoH < 24
        xline(ax, photoH, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    xlim(ax, [0 24]);
    script8_style_axes_(ax, theme);
end

function yLim = script8_activity_zt_ymax_(Act, clusterID, facets)
    yLim = [-0.5 1.5];
    if isempty(Act)
        return;
    end
    clusterIDs = string(clusterID);
    clusterIDs = clusterIDs(strlength(clusterIDs) > 0);
    if isempty(clusterIDs)
        return;
    end
    needed = {'ClusterID', 'Photoperiod_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames))
        return;
    end
    A = Act(ismember(string(Act.ClusterID), clusterIDs) & ismember(Act.Photoperiod_h, facets), :);
    if isempty(A), return; end
    y = double(A.Activity_zscored);
    y = y(isfinite(y));
    if isempty(y), return; end
    pad = 0.12 * max(range(y), 0.5);
    yLim = [min(y) - pad, max(y) + pad];
end

function hasData = script8_plot_zt_activity_facet_(ax, Act, clusterID, photo, pal, theme, yLim)
    if nargin < 7, yLim = []; end
    hasData = false;
    if isempty(Act) || strlength(string(clusterID)) == 0
        return;
    end
    needed = {'ClusterID', 'Photoperiod_h', 'SignalID', 'ZTBinCenter_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames))
        return;
    end
    A = Act(string(Act.ClusterID) == string(clusterID) & Act.Photoperiod_h == photo, :);
    if isempty(A)
        return;
    end
    hold(ax, 'on');
    set(ax, 'Color', 'w');
    if ~ismember('File', A.Properties.VariableNames)
        A.File = repmat("", height(A), 1);
    end
    mouseKey = string(A.File) + "|" + string(A.SignalID);
    sigs = unique(mouseKey, 'stable');
    mouseZT = {};
    mouseY = {};
    for s = 1:numel(sigs)
        As = A(mouseKey == sigs(s), :);
        As = sortrows(As, 'ZTBinCenter_h');
        zt = double(As.ZTBinCenter_h);
        y = double(As.Activity_zscored);
        keep = zt >= 0 & zt <= 24 & isfinite(y);
        if ~any(keep), continue; end
        plot(ax, zt(keep), y(keep), '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 0.85, 'HandleVisibility', 'off');
        mouseZT{end+1,1} = zt(keep); %#ok<AGROW>
        mouseY{end+1,1} = y(keep); %#ok<AGROW>
    end
    if ~isempty(mouseZT)
        hasData = true;
        edges = 0:0.5:24;
        centers = edges(1:end-1) + diff(edges)/2;
        meanY = nan(size(centers));
        for i = 1:numel(centers)
            vals = nan(numel(mouseZT), 1);
            for m = 1:numel(mouseZT)
                idx = mouseZT{m} >= edges(i) & mouseZT{m} < edges(i+1);
                if i == numel(centers)
                    idx = mouseZT{m} >= edges(i) & mouseZT{m} <= edges(i+1);
                end
                if any(idx)
                    vals(m) = mean(mouseY{m}(idx), 'omitnan');
                end
            end
            meanY(i) = mean(vals, 'omitnan');
        end
        plot(ax, centers, meanY, '-', 'Color', pal.base(1, :), 'LineWidth', 2.6);
    end
    if ~isempty(yLim)
        ylim(ax, yLim);
        script8_shade_zt_ld_(ax, double(photo), pal, yLim);
    else
        yl = ylim(ax);
        if diff(yl) < 1e-6
            yl = [-0.5 1.5];
            ylim(ax, yl);
        end
        script8_shade_zt_ld_(ax, double(photo), pal, yl);
    end
    xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    photoH = double(photo);
    if photoH > 0 && photoH < 24
        xline(ax, photoH, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
    xlim(ax, [0 24]);
    script8_style_axes_(ax, theme);
end

function script8_shade_zt_ld_(ax, photoH, pal, yl)
    if nargin < 4 || isempty(yl)
        yl = ylim(ax);
    end
    if yl(2) <= yl(1)
        yl = [-1 1];
    end
    photoH = double(photoH);
    if photoH >= 24 || photoH <= 0
        return;
    end
    h = patch(ax, [photoH 24 24 photoH], [yl(1) yl(1) yl(2) yl(2)], pal.cr, ...
        'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if ~isempty(h)
        uistack(h, 'bottom');
    end
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

function [manifest, standPaths] = script8_build_cluster_identity_cards_(data, outDirs, theme, manifest)
    standDir = fullfile(outDirs.standalone, 'Supplementary');
    standPaths = {};
    pal = theme.palette;
    identityPhotos = [12 22];
    nBands = min(2, numel(data.primaryUR));
    if nBands == 0
        return;
    end
    bands = data.primaryUR(1:nBands);
    clusters = strings(nBands, 1);
    for bi = 1:nBands
        clusters(bi) = script8_primary_cluster_(data.clusterSummary, bands(bi));
    end
    yMaxCoh = script8_coherence_zt_ymax_(data.phase24, clusters, identityPhotos, pal);
    yMaxAct = script8_activity_zt_ymax_(data.activityZT, clusters, identityPhotos);
    filterNotes = strings(0, 1);

    figAll = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1800 420 * nBands]);
    notesAll = script8_fill_identity_panels_(figAll, nBands, bands, clusters, identityPhotos, ...
        data, pal, theme, yMaxCoh, yMaxAct, 'ABCDEFGH', 0);
    filterNotes = unique([filterNotes; notesAll], 'stable');
    sgtitle(figAll, 'Primary-cluster identity cards (24h phase coherence + activity, L12 and L22)', ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
    baseAll = fullfile(standDir, 'Supp_ClusterIdentity_UR13_UR36');
    outAll = script8_export_figure_(figAll, baseAll, theme, {theme.ext, '.pdf'});
    close(figAll);
    standPaths{end + 1} = outAll{1}; %#ok<AGROW>
    wideCopy = fullfile(outDirs.compositeWide, ['Supp_ClusterIdentity_UR13_UR36' theme.ext]);
    copyfile(outAll{1}, wideCopy, 'f');
    cohFilterNote = 'Cluster-resolved PhaseCoherence_24h (ClusterID filter).';
    if any(filterNotes == "band-pooled (no ClusterID column)")
        cohFilterNote = 'Phase coherence band-pooled (PhaseCoherence_24h lacks ClusterID).';
    end
    manifest = script8_manifest_add_(manifest, 'Supp', 'ClusterIdentity', outAll{1}, 'standalone', ...
        ['Primary-cluster identity cards: 24h phase coherence and ZT activity at L12 and L22 for UR 1–3 and UR 3–6 primary clusters. ' ...
        'Supersedes band-pooled Supp_24h_L12_vs_L24.'], ...
        'Script 7', 'Tol band + LD shade', cohFilterNote);

    bandStems = ["UR13", "UR36"];
    for bi = 1:nBands
        figBand = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1800 420]);
        notesBand = script8_fill_identity_panels_(figBand, 1, bands(bi), clusters(bi), identityPhotos, ...
            data, pal, theme, yMaxCoh, yMaxAct, 'ABCD', 0);
        filterNotes = unique([filterNotes; notesBand], 'stable');
        sgtitle(figBand, sprintf('Cluster identity — %s (%s)', ...
            script8_band_display_(bands(bi), 'plain'), script8_cluster_short_label_(clusters(bi), data.clusterSummary)), ...
            'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
        baseBand = fullfile(standDir, ['Supp_ClusterIdentity_' bandStems(bi)]);
        outBand = script8_export_figure_(figBand, baseBand, theme, {theme.ext, '.pdf'});
        close(figBand);
        standPaths{end + 1} = outBand{1}; %#ok<AGROW>
        manifest = script8_manifest_add_(manifest, 'Supp', ['ClusterIdentity_' char(bandStems(bi))], outBand{1}, 'standalone', ...
            sprintf('Primary-cluster identity card for %s: L12/L22 coherence and activity.', script8_band_display_(bands(bi), 'plain')), ...
            'Script 7', 'Tol band + LD shade', cohFilterNote);
    end
end

function notes = script8_fill_identity_panels_(fig, nRows, rowBands, rowClusters, identityPhotos, data, pal, theme, yMaxCoh, yMaxAct, panelLetters, letterOffset)
    notes = strings(0, 1);
    rowBands = string(rowBands);
    rowClusters = string(rowClusters);
    for ri = 1:nRows
        bn = rowBands(ri);
        cid = rowClusters(ri);
        bandCol = script8_band_colour_(pal, bn);
        cidShort = script8_cluster_short_label_(cid, data.clusterSummary);
        for pj = 1:numel(identityPhotos)
            photo = identityPhotos(pj);
            ppTxt = char(script8_pp_label_(photo));
            tileCoh = (ri - 1) * 4 + (pj - 1) * 2 + 1;
            tileAct = tileCoh + 1;
            axC = subplot(nRows, 4, tileCoh, 'Parent', fig);
            [hasCoh, noteCoh] = script8_plot_zt_coherence_facet_(axC, data.phase24, cid, photo, bandCol, pal, theme, yMaxCoh, bn);
            if strlength(noteCoh) > 0, notes(end + 1, 1) = noteCoh; end %#ok<AGROW>
            if ~hasCoh
                text(axC, 0.5, 0.5, sprintf('No coherence (%s)', ppTxt), 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'FontName', theme.fontName);
            end
            title(axC, sprintf('%s | %s | %s coherence', script8_band_display_(bn, 'plain'), cidShort, ppTxt), ...
                'FontWeight', 'bold', 'Interpreter', 'none');
            if pj == 1
                ylabel(axC, 'Phase coherence R', 'FontWeight', 'bold');
            end
            if ri == nRows
                xlabel(axC, 'ZT (h)', 'FontWeight', 'bold');
            end
            li = letterOffset + (ri - 1) * 4 + (pj - 1) * 2 + 1;
            if li <= numel(panelLetters)
                script8_panel_label_(axC, panelLetters(li), theme);
            end

            axA = subplot(nRows, 4, tileAct, 'Parent', fig);
            hasAct = script8_plot_zt_activity_facet_(axA, data.activityZT, cid, photo, pal, theme, yMaxAct);
            if ~hasAct
                text(axA, 0.5, 0.5, sprintf('No activity (%s)', ppTxt), 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'FontName', theme.fontName);
            end
            title(axA, sprintf('%s | %s | %s activity', script8_band_display_(bn, 'plain'), cidShort, ppTxt), ...
                'FontWeight', 'bold', 'Interpreter', 'none');
            if pj == 1
                ylabel(axA, 'Activity (z-scored)', 'FontWeight', 'bold');
            end
            if ri == nRows
                xlabel(axA, 'ZT (h)', 'FontWeight', 'bold');
            end
            li = li + 1;
            if li <= numel(panelLetters)
                script8_panel_label_(axA, panelLetters(li), theme);
            end
        end
    end
end

%% Supplementary — sex + cluster identity + non-primary activity + legacy + coexpression
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

    % Primary-cluster identity cards (replaces band-pooled Supp_24h_L12_vs_L24)
    [manifest, idPaths] = script8_build_cluster_identity_cards_(data, outDirs, theme, manifest);
    standPaths = [standPaths, idPaths]; %#ok<AGROW>

    % Per-cluster 24h activity (non-primary UR 1–3 clusters)
    facets = pal.coherenceFacets;
    CS = data.clusterSummary;
    primaryCluster = script8_primary_cluster_(CS, data.primaryUR(1));
    if ~isempty(CS) && ismember('ClusterID', CS.Properties.VariableNames)
        for ci = 1:height(CS)
            cid = string(CS.ClusterID(ci));
            if strlength(cid) == 0 || cid == primaryCluster
                continue;
            end
            yMaxAct = script8_activity_zt_ymax_(data.activityZT, cid, facets);
            figAct = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1400 780]);
            for fi = 1:numel(facets)
                axS = subplot(2, 3, fi); hold(axS, 'on');
                set(axS, 'Color', 'w');
                hasData = script8_plot_zt_activity_facet_(axS, data.activityZT, cid, facets(fi), pal, theme, yMaxAct);
                title(axS, char(script8_pp_label_(facets(fi))), 'FontWeight', 'bold');
                if ~hasData
                    text(axS, 0.5, 0.5, sprintf('No activity for %s', char(script8_pp_label_(facets(fi)))), ...
                        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontName', theme.fontName);
                end
                if fi == 1 || fi == 4
                    ylabel(axS, 'Activity (z-scored)', 'FontWeight', 'bold');
                end
                if fi >= 4
                    xlabel(axS, 'ZT (h)', 'FontWeight', 'bold');
                end
                set(axS, 'XLim', [0 24], 'Color', 'w');
                script8_style_axes_(axS, theme);
            end
            cidSafe = regexprep(char(cid), '[^A-Za-z0-9_]+', '_');
            sgtitle(figAct, sprintf('24h activity — cluster %s', char(cid)), ...
                'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
            outAct = script8_export_figure_(figAct, fullfile(standDir, ['Supp_Activity_L12_L22_' cidSafe]), theme, {theme.ext, '.pdf'});
            close(figAct);
            standPaths{end + 1} = outAct{1}; %#ok<AGROW>
            manifest = script8_manifest_add_(manifest, 'Supp', ['Activity_' cidSafe], outAct{1}, 'standalone', ...
                sprintf('24h ZT activity L12–L22 for cluster %s (non-primary).', char(cid)), ...
                'Script 7', 'Tol L12/L22', 'Excluded primary UR 1–3 cluster (on Fig06).');
        end
    end

    % Legacy dual-band coherence / ridge-at-transition (after identity cards)
    figLegacy = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1200 780]);
    heroPP = [12 16 22];
    legacyYMax = script8_coherence_ymax_(data.binnedCoherence, heroPP, "UR_1_3", pal);
    legacyYMax = max(legacyYMax, script8_coherence_ymax_(data.binnedCoherence, heroPP, "UR_3_6", pal));
    for fi = 1:numel(heroPP)
        axB = subplot(2, 3, fi); hold(axB, 'on'); set(axB, 'Color', 'w');
        script8_plot_coherence_facet_(axB, data.binnedCoherence, heroPP(fi), "UR_1_3", pal, theme, false);
        title(axB, [char(script8_pp_label_(heroPP(fi))) ' UR 1–3'], 'FontWeight', 'bold');
        set(axB, 'YLim', [0 legacyYMax], 'XLim', pal.coherenceXlim);
        script8_style_axes_(axB, theme);
        axC = subplot(2, 3, fi + 3); hold(axC, 'on'); set(axC, 'Color', 'w');
        script8_plot_coherence_facet_(axC, data.binnedCoherence, heroPP(fi), "UR_3_6", pal, theme, false);
        title(axC, [char(script8_pp_label_(heroPP(fi))) ' UR 3–6'], 'FontWeight', 'bold');
        set(axC, 'YLim', [0 legacyYMax], 'XLim', pal.coherenceXlim);
        script8_style_axes_(axC, theme);
    end
    outL = script8_export_figure_(figLegacy, fullfile(standDir, 'Supp_Legacy_DualBand_Coherence'), theme, {theme.ext, '.pdf'});
    close(figLegacy);
    standPaths{end + 1} = outL{1}; %#ok<AGROW>
    manifest = script8_manifest_add_(manifest, 'Supp', 'Legacy_coherence', outL{1}, 'standalone', ...
        'Former main dual-band coherence facets (now Supplementary).', 'Scripts 5/7', 'Tol', '');

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
        'Former main ridge-at-transition panels (now Supplementary).', 'Scripts 5/7', 'Tol', '');

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
        'L12 (compare)', '#4477AA', 'Fig06 / supp 24h activity + cluster identity L12';
        'L22 / L24', '#AA3377', 'Fig06 / supp 24h activity + cluster identity L22';
        'Cluster identity coherence', 'Band colour', 'Supp_ClusterIdentity_* coherence traces';
        'LD shading (ZT)', '#BBBBBB 8%', '24h activity + cluster identity ZT panels';
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
    set(ax, 'Box', 'off', 'TickDir', theme.palette.tickDir, 'FontName', theme.fontName, ...
        'LineWidth', theme.palette.axesLineWidth, 'FontSize', 11);
    if isprop(ax, 'XLabel') && ~isempty(ax.XLabel.String)
        ax.XLabel.FontWeight = 'bold';
        ax.XLabel.FontSize = 12;
    end
    if isprop(ax, 'YLabel') && ~isempty(ax.YLabel.String)
        ax.YLabel.FontWeight = 'bold';
        ax.YLabel.FontSize = 12;
    end
end

function script8_panel_label_(ax, labelChar, theme)
    text(ax, 0.02, 0.98, labelChar, 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 16, 'FontName', theme.fontName, 'VerticalAlignment', 'top');
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
