function out = extended_script10_manuscript_figures_run(cohortRoot, cfg)
%EXTENDED_SCRIPT10_MANUSCRIPT_FIGURES_RUN Locked manuscript figure set (Figs 1–7).
%
%   Fig01 Circadian characteristics — placeholder (collaborator external)
%   Fig02 Multiscale RAW + HSub residual Female|Male + CarryForward retention E
%   Fig03 Co-expression: A absolute power; B descriptive Delta; C LME forest
%   Fig04 UR_1_3 Z-scored activity L12–L22 (C01 / ClusterRank 1)
%   Fig05 UR_1_3 24h ZT coherence grid only (main composite = panel A)
%   Fig06 UR_3_6 activity twin of Fig04 (C01)
%   Fig07 UR_3_6 coherence twin of Fig05 (C01; A-only main composite)
%
%   Supplement (not main storyboard / not renumbered mains):
%     Standalone/Supp_Transitions/ — LD-only paired Pre/Post phase-coherence R bars
%       stems: Supp_Fig05_LD_PrePost_R, Supp_Fig07_LD_PrePost_R (UR_1_3 / UR_3_6)
%       BH-FDR stars from Resync_PrimaryStats_BH_FDR (DeltaR>0); cite Script 5 tables
%     Standalone/Supp/ — UR_3_6 C02 activity + coherence twins (cfg manuscriptClusters)
%       stems: Supp_Activity_UR36_C02, Supp_Coherence_UR36_C02
%
%   Cluster lock: cfg.script10.manuscriptClusters (ClusterRank; optional period check).
%   Sex after Fig02 stays Script 9 — not in main composites.
%   Package_Manuscript/ holds Figs 3–7 + Tables (full supplemental export) + Figures/Supp/ (excludes Fig01 & Fig02).
%
%   Does NOT invent new ridge/LME maths — reglues Scripts 5–9 outputs only.

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    if ~isfield(cfg, 'plotMode') || isempty(cfg.plotMode)
        cfg.plotMode = 'development';
    end
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);
    theme = script10_theme_(cfg);

    if nargin < 1 || isempty(cohortRoot)
        cohortRoot = uigetdir(pwd, 'Select cohort results folder (e.g. C57_LP)');
        if isequal(cohortRoot, 0)
            error('extended_script10_manuscript_figures_run:NoRoot', 'No cohort folder selected.');
        end
    end

    paths = extended_script8_resolve_paths(cohortRoot);
    if ~isempty(paths.missing)
        error('extended_script10_manuscript_figures_run:MissingInputs', ...
            'Missing required inputs: %s', strjoin(paths.missing, ', '));
    end

    [outDirs, modeLabel] = script10_make_output_dirs_(cohortRoot, cfg.plotMode);
    logPath = fullfile(outDirs.logs, sprintf('Script10_ManuscriptFigures_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() script10_fclose_(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 10: Manuscript figures ===\n');
    fprintf('Cohort:  %s\n', paths.cohortRoot);
    fprintf('Output:  %s\n', outDirs.root);
    fprintf('Mode:    %s → folder %s | dpi=%g | ext=%s\n', ...
        char(string(cfg.plotMode)), char(modeLabel), theme.dpi, theme.ext);
    script10_log_(LOG, 'Cohort: %s', paths.cohortRoot);
    script10_log_(LOG, 'plotMode=%s | folder=%s | dpi=%g | ext=%s', ...
        char(string(cfg.plotMode)), char(modeLabel), theme.dpi, theme.ext);

    data = script10_load_data_(paths, cfg);
    manifest = script10_manifest_init_();
    compositeWide = strings(0, 1);
    packageFigs = strings(0, 1);  % Fig03–Fig07 only (main composites)
    packageSupp = strings(0, 1); % Standalone/Supp* paths for Package Figures/Supp/

    %% Fig01 — Circadian placeholder (collaborator external)
    [manifest, wide1, ~] = script10_build_fig01_placeholder_(outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide1); %#ok<AGROW>
    script10_log_(LOG, 'Fig 01 placeholder complete');

    %% Fig02 — RAW + HSub Female|Male + retention E
    [manifest, wide2, tall2, ~] = script10_build_fig02_(data, paths, outDirs, theme, manifest);
    compositeWide(end + 1, 1) = string(wide2); %#ok<AGROW>
    script10_export_tall_copy_(tall2, wide2, outDirs.compositeTall);
    script10_log_(LOG, 'Fig 02 complete');

    %% Fig03 — Co-expression A abs power | B descriptive Delta | C LME forest
    [manifest, wide3, tall3] = script10_build_fig03_(data, outDirs, theme, cfg, manifest);
    compositeWide(end + 1, 1) = string(wide3); %#ok<AGROW>
    packageFigs(end + 1, 1) = string(wide3); %#ok<AGROW>
    script10_export_tall_copy_(tall3, wide3, outDirs.compositeTall);
    script10_log_(LOG, 'Fig 03 complete');

    %% Fig04 — UR_1_3 activity L12–L22 (C01 only)
    band13 = data.primaryUR(1);
    mainRank13 = script10_manuscript_main_rank_(cfg, band13);
    cid13 = script10_resolve_manuscript_cluster_(data.clusterSummary, band13, mainRank13, cfg, LOG);
    [manifest, wide4, tall4] = script10_build_activity_figure_(data, outDirs, theme, cfg, manifest, ...
        'Fig04', band13, 'Fig04_Activity_UR13', cid13, true, '');
    compositeWide(end + 1, 1) = string(wide4); %#ok<AGROW>
    packageFigs(end + 1, 1) = string(wide4); %#ok<AGROW>
    script10_export_tall_copy_(tall4, wide4, outDirs.compositeTall);
    script10_log_(LOG, 'Fig 04 complete (cluster %s)', char(cid13));

    %% Fig05 — UR_1_3 coherence A-only main + LD Pre/Post R supp
    [manifest, wide5, tall5, supp5] = script10_build_coherence_figure_(data, outDirs, theme, cfg, manifest, ...
        'Fig05', band13, 'Fig05_Coherence_UR13', cid13, true, true, '');
    compositeWide(end + 1, 1) = string(wide5); %#ok<AGROW>
    packageFigs(end + 1, 1) = string(wide5); %#ok<AGROW>
    packageSupp = [packageSupp; supp5(:)]; %#ok<AGROW>
    script10_export_tall_copy_(tall5, wide5, outDirs.compositeTall);
    script10_log_(LOG, 'Fig 05 complete (A-only main; LD Pre/Post R → Supp_Transitions)');

    %% Fig06 — UR_3_6 activity twin (C01)
    band36 = data.primaryUR(min(2, numel(data.primaryUR)));
    mainRank36 = script10_manuscript_main_rank_(cfg, band36);
    cid36 = script10_resolve_manuscript_cluster_(data.clusterSummary, band36, mainRank36, cfg, LOG);
    [manifest, wide6, tall6] = script10_build_activity_figure_(data, outDirs, theme, cfg, manifest, ...
        'Fig06', band36, 'Fig06_Activity_UR36', cid36, true, '');
    compositeWide(end + 1, 1) = string(wide6); %#ok<AGROW>
    packageFigs(end + 1, 1) = string(wide6); %#ok<AGROW>
    script10_export_tall_copy_(tall6, wide6, outDirs.compositeTall);
    script10_log_(LOG, 'Fig 06 complete (cluster %s)', char(cid36));

    %% Fig07 — UR_3_6 coherence twin (C01; A-only + LD Pre/Post R supp)
    [manifest, wide7, tall7, supp7] = script10_build_coherence_figure_(data, outDirs, theme, cfg, manifest, ...
        'Fig07', band36, 'Fig07_Coherence_UR36', cid36, true, true, '');
    compositeWide(end + 1, 1) = string(wide7); %#ok<AGROW>
    packageFigs(end + 1, 1) = string(wide7); %#ok<AGROW>
    packageSupp = [packageSupp; supp7(:)]; %#ok<AGROW>
    script10_export_tall_copy_(tall7, wide7, outDirs.compositeTall);
    script10_log_(LOG, 'Fig 07 complete (A-only main; LD Pre/Post R → Supp_Transitions)');

    %% Supp extras — UR_3_6 C02 (rank 2) activity + coherence (not main Figs / storyboard)
    suppRanks36 = script10_manuscript_supp_ranks_(cfg, band36);
    for ri = 1:numel(suppRanks36)
        rnk = double(suppRanks36(ri));
        cidSupp = script10_resolve_manuscript_cluster_(data.clusterSummary, band36, rnk, cfg, LOG);
        if strlength(string(cidSupp)) == 0
            script10_log_(LOG, 'Supp skip: no ClusterRank %g for %s', rnk, char(string(band36)));
            continue;
        end
        short = script10_cluster_short_label_(cidSupp, data.clusterSummary);
        actStem = sprintf('Supp_Activity_UR36_%s', short);
        cohStem = sprintf('Supp_Coherence_UR36_%s', short);
        [manifest, ~, ~] = script10_build_activity_figure_(data, outDirs, theme, cfg, manifest, ...
            'Supp', band36, actStem, cidSupp, false, 'Supp');
        packageSupp(end + 1, 1) = string(fullfile(outDirs.standalone, 'Supp', actStem)); %#ok<AGROW>
        [manifest, ~, ~, ~] = script10_build_coherence_figure_(data, outDirs, theme, cfg, manifest, ...
            'Supp', band36, cohStem, cidSupp, false, false, 'Supp');
        packageSupp(end + 1, 1) = string(fullfile(outDirs.standalone, 'Supp', cohStem)); %#ok<AGROW>
        script10_log_(LOG, 'Supp UR_3_6 %s activity+coherence exported (not main)', short);
    end

    script10_write_manifest_(manifest, outDirs.manifest);
    script10_write_figure_manifest_md_(outDirs.root, manifest);
    script10_write_figure_manifest_csv_(outDirs.root, manifest);
    script10_save_storyboard_(compositeWide, outDirs.storyboard);

    %% Package handoff (Figs 3–7 main + Supp folder + tables; exclude Fig01 & Fig02)
    pkgRoot = script10_build_package_(paths, outDirs, packageFigs, packageSupp, theme, cfg, LOG);
    script10_log_(LOG, 'Package_Manuscript: %s', pkgRoot);

    out = struct();
    out.cohortRoot = paths.cohortRoot;
    out.outRoot = outDirs.root;
    out.packageRoot = pkgRoot;
    out.manifest = outDirs.manifest;
    out.storyboard = outDirs.storyboard;
    out.nManifestRows = height(manifest);
    out.logPath = logPath;

    fprintf('\nExtended Script 10 complete.\n');
    fprintf('  Composites: %s\n', outDirs.compositeWide);
    fprintf('  Package:    %s\n', pkgRoot);
    fprintf('  Manifest:   %s\n', outDirs.manifest);
    fprintf(['\nNote: Sex-stratified and other supplements remain Script 9.\n' ...
        '  Run:  run_extended_script9_supplementary_figures\n']);
end

%% ========================================================================
function data = script10_load_data_(paths, cfg)
    data = struct();
    data.pairSummary = script10_read_sheet_(paths.lmeDescriptiveXlsx, 'CR_UR_Pairs_Summary');
    data.absSummary = script10_read_sheet_(paths.lmeDescriptiveXlsx, 'AbsolutePower_Summary');
    data.sexBalance = script10_read_sheet_(paths.lmeDescriptiveXlsx, 'Sex_Balance');
    data.lmeCoefDelta = script10_read_sheet_(paths.lmeInferenceXlsx, 'LME_Coef_Delta_BH_FDR');
    data.lmeCoefPower = script10_read_sheet_(paths.lmeInferenceXlsx, 'LME_Coef_Power_BH_FDR');
    data.lmeFdr = script10_read_sheet_(paths.lmeInferenceXlsx, 'LME_Inference_BH_FDR');
    data.binnedCoherence = script10_read_sheet_(paths.resyncXlsx, 'BinnedCoherence');
    data.candidateDeltaR = script10_read_sheet_(paths.resyncXlsx, 'CandidateDeltaR');
    data.deltaRSummary = script10_read_sheet_(paths.resyncXlsx, 'DeltaR_Summary');
    data.resyncPrimaryFdr = script10_read_sheet_(paths.resyncXlsx, 'Resync_PrimaryStats_BH_FDR');
    data.clusterSummary = script10_read_sheet_(paths.profilesXlsx, 'ClusterSummary');
    data.activityZT = script10_read_sheet_(paths.profilesXlsx, 'ActivityComponent_24h');
    data.phase24 = script10_read_sheet_(paths.profilesXlsx, 'PhaseCoherence_24h');
    data.retention = script10_read_sheet_(paths.gateXlsx, 'Retention_ByBand');
    data.primaryUR = string(cfg.bands.primaryUR);
    data.ppOrder = script10_pp_order_(data.pairSummary);
    data.nMice = script10_total_mice_(data.sexBalance);
end

%% Fig01 — Circadian characteristics placeholder
function [manifest, widePath, standPaths] = script10_build_fig01_placeholder_(outDirs, theme, manifest)
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig01_CircadianPlaceholder');
    extended_period_gate_ensure_dir(standDir);

    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 1400 520], 'Visible', 'off');
    ax = axes(fig, 'Position', [0.06 0.12 0.88 0.76]); hold(ax, 'on');
    axis(ax, [0 10 0 1]); axis(ax, 'off');

    rectangle(ax, 'Position', [0.8 0.25 8.4 0.55], 'Curvature', 0.06, ...
        'FaceColor', [0.97 0.97 0.97], 'EdgeColor', pal.base(1, :), 'LineWidth', 2.0, 'LineStyle', '--');
    text(ax, 5, 0.68, 'FIGURE 1 — PLACEHOLDER', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 20, 'Color', pal.base(1, :));
    text(ax, 5, 0.50, 'Circadian characteristics (collaborator external figure)', ...
        'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontSize', 14, 'Interpreter', 'none');
    text(ax, 5, 0.36, 'Drop final circadian panel assets here before manuscript lock.', ...
        'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontSize', 11, ...
        'Color', [0.35 0.35 0.35], 'Interpreter', 'none');

    base = fullfile(standDir, 'Fig01_Circadian_Placeholder');
    standPaths = script10_export_figure_(fig, base, theme, {theme.ext, '.pdf'});
    widePath = fullfile(outDirs.compositeWide, ['Fig01_Circadian_Placeholder' theme.ext]);
    copyfile(standPaths{1}, widePath, 'f');
    close(fig);

    readmePath = fullfile(standDir, 'README_Fig01_SLOT.md');
    fid = fopen(readmePath, 'w');
    if fid > 0
        fprintf(fid, '# Fig 1 — Circadian characteristics (slot)\n\n');
        fprintf(fid, 'This panel is **not generated by Extended Scripts 4–10**.\n\n');
        fprintf(fid, 'Expected content: collaborator-provided circadian characterisation figure.\n\n');
        fprintf(fid, 'Replace `Fig01_Circadian_Placeholder` assets in Standalone and CompositePanels/Wide_16x9\n');
        fprintf(fid, 'before manuscript submission. Keep filename stem `Fig01_*` for storyboard order.\n\n');
        fprintf(fid, '**Excluded from Package_Manuscript/** (collaborator-owned; not part of Extended handoff zip).\n');
        fclose(fid);
    end

    manifest = script10_manifest_add_(manifest, 'Fig01', 'Placeholder', widePath, '16:9', ...
        'PLACEHOLDER: Circadian characteristics — collaborator external figure (not Extended pipeline).', ...
        'External', 'Schematic', 'Exploratory slot only; replace before lock. Excluded from Package_Manuscript.');
end

%% Fig02 — RAW A|B, HSub C|D, retention E
function [manifest, widePath, tallPath, standPaths] = script10_build_fig02_(data, paths, outDirs, theme, manifest)
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
            error('extended_script10_manuscript_figures_run:MissingScalogram', ...
                'Missing scalogram image: %s', panels{i, 1});
        end
        if panels{i, 4}
            img = script10_crop_hsub_residual_(panels{i, 1});
        else
            img = imread(panels{i, 1});
        end
        fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'pixels', ...
            'Position', [50 50 max(size(img, 2), 400) max(size(img, 1), 300)]);
        ax = axes(fig); imshow(img, 'Parent', ax); axis(ax, 'off');
        title(ax, panels{i, 3}, 'FontName', theme.fontName, 'FontWeight', 'bold', 'Interpreter', 'none');
        script10_panel_label_(ax, panels{i, 2}, theme);
        base = fullfile(standDir, sprintf('Fig02_%s_%s', panels{i, 2}, ...
            regexprep(panels{i, 3}, '[^A-Za-z0-9]+', '_')));
        standPaths{i} = script10_export_figure_(fig, base, theme, {theme.ext});
        close(fig);
        note = 'RAW circadian context; Female|Male not pooled.';
        if panels{i, 4}, note = 'HSub residual half-crop; Female|Male not pooled.'; end
        manifest = script10_manifest_add_(manifest, 'Fig02', panels{i, 2}, standPaths{i}{1}, 'standalone', ...
            panels{i, 3}, 'Script 3', 'Jet scalogram', note);
    end

    figE = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 720 420]);
    axE = axes(figE); hold(axE, 'on');
    R = data.retention;
    if ~isempty(R)
        bn = string(R.BandName);
        retCol = R.Properties.VariableNames{contains(R.Properties.VariableNames, 'RetentionPct_Elig', 'IgnoreCase', true)};
        y = double(R.(retCol));
        for i = 1:numel(y)
            bar(axE, i, y(i), 0.7, 'FaceColor', script10_band_colour_(pal, bn(i)), 'EdgeColor', 'k');
        end
        set(axE, 'XTick', 1:numel(bn), 'XTickLabel', arrayfun(@(b) script10_band_display_(b, 'tex'), bn, 'UniformOutput', false));
        axE.TickLabelInterpreter = 'tex';
        ylim(axE, [0 105]);
    end
    ylabel(axE, 'Validated retention (% eligible Raw)', 'FontWeight', 'bold');
    title(axE, 'CarryForward retention by band', 'FontWeight', 'bold');
    script10_style_axes_(axE, theme);
    script10_panel_label_(axE, 'E', theme);
    standPaths{5} = script10_export_figure_(figE, fullfile(standDir, 'Fig02_E_Retention'), theme, {theme.ext, '.pdf'});
    close(figE);
    manifest = script10_manifest_add_(manifest, 'Fig02', 'E', standPaths{5}{1}, 'standalone', ...
        'CarryForward retention by band (QC)', 'Script 4', 'Tol bands', 'Optional QC panel.');

    [widePath, tallPath] = script10_composite_fig02_(standPaths, outDirs, theme);
    manifest = script10_manifest_add_(manifest, 'Fig02', 'Composite', widePath, '16:9', ...
        'Multiscale RAW + HSub residual scalograms (Female|Male) with CarryForward retention QC.', ...
        'Scripts 3/4', 'Jet + Tol', 'Excluded from Package_Manuscript.');
end

function [widePath, tallPath] = script10_composite_fig02_(standPaths, outDirs, theme)
    imgs = cell(5, 1);
    for i = 1:5
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script10_trim_image_whitespace_(imread(p));
    end
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 1100]);
    tl = tiledlayout(figW, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1); imshow(imgs{1}); axis off;
    nexttile(tl, 2); imshow(imgs{2}); axis off;
    nexttile(tl, 3); imshow(imgs{3}); axis off;
    nexttile(tl, 4); imshow(imgs{4}); axis off;
    nexttile(tl, 5, [1 2]); imshow(imgs{5}); axis off;
    widePath = fullfile(outDirs.compositeWide, ['Fig02_Scalograms_Retention' theme.ext]);
    script10_exportgraphics_(figW, widePath, theme.dpi);
    close(figW);

    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1400]);
    tl2 = tiledlayout(figT, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl2, 1); imshow(imgs{1}); axis off;
    nexttile(tl2, 2); imshow(imgs{2}); axis off;
    nexttile(tl2, 3); imshow(imgs{3}); axis off;
    nexttile(tl2, 4); imshow(imgs{4}); axis off;
    nexttile(tl2, 5, [1 2]); imshow(imgs{5}); axis off;
    tallPath = fullfile(outDirs.compositeTall, ['Fig02_Scalograms_Retention' theme.ext]);
    script10_exportgraphics_(figT, tallPath, theme.dpi);
    close(figT);
end

%% Fig03 — A absolute power | B descriptive Delta | C LME forest
function [manifest, widePath, tallPath] = script10_build_fig03_(data, outDirs, theme, cfg, manifest) %#ok<INUSD>
    pal = theme.palette;
    standDir = fullfile(outDirs.standalone, 'Fig03_Coexpression');
    extended_period_gate_ensure_dir(standDir);
    standPaths = cell(3, 1);
    pp = data.ppOrder;
    ppLabels = arrayfun(@(x) char(script10_pp_label_(x)), pp, 'UniformOutput', false);

    %% A — Absolute power CR + primary URs
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axA = axes(figA); hold(axA, 'on');
    cr = data.absSummary(string(data.absSummary.BandName) == "CR_20_28" & string(data.absSummary.Phase) == "All", :);
    cr = sortrows(cr, 'Photoperiod_h');
    if ~isempty(cr)
        crErr = cr.SD_Log10 ./ sqrt(max(data.nMice, 1));
        errorbar(axA, 1:height(cr), cr.Mean_Log10, crErr, '-s', 'Color', pal.cr, 'LineWidth', 1.5, ...
            'MarkerFaceColor', pal.cr, 'CapSize', 6);
        script10_direct_line_label_(axA, height(cr), cr.Mean_Log10(end), script10_band_display_('CR_20_28', 'tex'), pal.cr, theme, true);
    end
    for b = 1:numel(data.primaryUR)
        bn = data.primaryUR(b);
        ur = data.absSummary(string(data.absSummary.BandName) == bn & string(data.absSummary.Phase) == "All", :);
        ur = sortrows(ur, 'Photoperiod_h');
        if isempty(ur), continue; end
        col = script10_band_colour_(pal, bn);
        urErr = ur.SD_Log10 ./ sqrt(max(data.nMice, 1));
        errorbar(axA, 1:height(ur), ur.Mean_Log10, urErr, '-o', 'Color', col, 'LineWidth', 2, ...
            'MarkerFaceColor', col, 'CapSize', 8);
        script10_direct_line_label_(axA, height(ur), ur.Mean_Log10(end), script10_band_display_(bn, 'tex'), col, theme, true);
    end
    set(axA, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axA, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axA, 'Mean log_{10} band power', 'FontWeight', 'bold');
    title(axA, 'Absolute power (CR, UR 1–3, UR 3–6)', 'FontWeight', 'bold');
    script10_style_axes_(axA, theme);
    script10_panel_label_(axA, 'A', theme);
    script10_add_n_annotation_(axA, data.nMice, theme);
    standPaths{1} = script10_export_figure_(figA, fullfile(standDir, 'Fig03_A_AbsPower'), theme, {theme.ext, '.pdf'});
    close(figA);

    %% B — Descriptive Delta = UR - CR (both primary bands)
    figB = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 620 480]);
    axB = axes(figB); hold(axB, 'on');
    for b = 1:numel(data.primaryUR)
        bn = data.primaryUR(b);
        sub = data.pairSummary(string(data.pairSummary.UR_Band) == bn & string(data.pairSummary.Phase) == "All", :);
        sub = sortrows(sub, 'Photoperiod_h');
        if isempty(sub), continue; end
        col = script10_band_colour_(pal, bn);
        y = sub.Mean_Delta_log10;
        x = 1:height(sub);
        err = sub.SD_Delta_log10 ./ sqrt(max(data.nMice, 1));
        errorbar(axB, x, y, err, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, 'CapSize', 8);
        script10_direct_line_label_(axB, x(end), y(end), script10_band_display_(bn, 'tex'), col, theme, true);
    end
    script10_yline_zero_(axB);
    set(axB, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
    xlabel(axB, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(axB, '\Delta(UR - CR) log_{10} power', 'FontWeight', 'bold');
    title(axB, 'Descriptive \Delta = UR − CR (primary bands)', 'FontWeight', 'bold');
    script10_style_axes_(axB, theme);
    script10_panel_label_(axB, 'B', theme);
    script10_add_n_annotation_(axB, data.nMice, theme);
    standPaths{2} = script10_export_figure_(figB, fullfile(standDir, 'Fig03_B_DeltaDescriptive'), theme, {theme.ext, '.pdf'});
    close(figB);

    %% C — LME forest (Photoperiod_h beta ± CI); confirmatory co-expression inference
    figC = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 720 480]);
    axC = axes(figC); hold(axC, 'on');
    script10_plot_lme_forest_(axC, data, theme);
    title(axC, 'LME Photoperiod_h effects (BH-FDR)', 'FontWeight', 'bold');
    script10_style_axes_(axC, theme);
    script10_panel_label_(axC, 'C', theme);
    standPaths{3} = script10_export_figure_(figC, fullfile(standDir, 'Fig03_C_LME_Forest'), theme, {theme.ext, '.pdf'});
    close(figC);

    captions = {
        'Absolute band power across photoperiod (descriptive; CR + primary URs).';
        'Descriptive \Delta=UR−CR for UR_1_3 and UR_3_6 (exploratory summary means ± SEM).';
        'Confirmatory LME Photoperiod_h \beta\pm95% CI from Script 6 (Delta primary; optional Power). Sex covariate OK; AdditiveOnly preferred.'};
    sources = {'Script 6', 'Script 6', 'Script 6'};
    for k = 1:3
        manifest = script10_manifest_add_(manifest, 'Fig03', char('A' + k - 1), standPaths{k}{1}, 'standalone', ...
            captions{k}, sources{k}, 'Tol', '');
    end
    [widePath, tallPath] = script10_composite_fig03_(standPaths, outDirs, theme);
    manifest = script10_manifest_add_(manifest, 'Fig03', 'Composite', widePath, '16:9', ...
        'CR–UR co-expression: absolute power (A), descriptive \Delta (B), LME forest (C). Sex → Script 9.', ...
        'Script 6', 'Tol', 'Fig3C = co-expression LME; do not conflate with Figs 5/7 transition FDR.');
end

function script10_plot_lme_forest_(ax, data, theme)
    pal = theme.palette;
    rows = script10_lme_forest_rows_(data);
    if isempty(rows)
        text(ax, 0.5, 0.5, 'No LME Photoperiod_h coefficients found', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        axis(ax, 'off');
        return;
    end

    n = height(rows);
    % Plot top-to-bottom: row 1 at highest y; YTick still ascending
    y = (n:-1:1)';
    hold(ax, 'on');
    xline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');

    for i = 1:n
        est = rows.Estimate(i);
        lo = rows.CI_lo(i);
        hi = rows.CI_hi(i);
        col = script10_band_colour_(pal, rows.BandName(i));
        if rows.MetricClass(i) == "Power"
            col = min(col * 0.75 + 0.15, 1);
        end
        plot(ax, [lo hi], [y(i) y(i)], '-', 'Color', col, 'LineWidth', 2.0);
        ms = 7;
        if rows.Significant_BH(i)
            plot(ax, est, y(i), 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', ms + 2, 'LineWidth', 1.2);
            text(ax, hi + 0.01 * max(1, abs(hi)), y(i), '*', 'Color', 'k', 'FontSize', 14, ...
                'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
        else
            plot(ax, est, y(i), 'o', 'Color', col, 'MarkerFaceColor', 'w', 'MarkerSize', ms, 'LineWidth', 1.4);
        end
    end

    labels = arrayfun(@(i) sprintf('%s  %s', char(rows.MetricClass(i)), ...
        script10_band_display_(rows.BandName(i), 'plain')), (1:n)', 'UniformOutput', false);
    % Ascending ticks with labels flipped to match top-to-bottom row order
    set(ax, 'YTick', 1:n, 'YTickLabel', flipud(labels), 'YLim', [0.4 n + 0.6], ...
        'TickLabelInterpreter', 'none');
    xlabel(ax, 'Photoperiod_h \beta (95% CI)', 'FontWeight', 'bold', 'Interpreter', 'tex');
    ylabel(ax, '');
    text(ax, 0.98, 0.02, '* Significant_BH (Script 6 co-expression LME)', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35], 'Interpreter', 'none');
end

function rows = script10_lme_forest_rows_(data)
%SCRIPT10_LME_FOREST_ROWS_ Prefer AdditiveOnly Photoperiod_h for primary UR Delta then Power.
    emptyVars = {'MetricClass','BandName','Estimate','SE','CI_lo','CI_hi','Significant_BH'};
    rows = table('Size', [0 7], ...
        'VariableTypes', {'string','string','double','double','double','double','logical'}, ...
        'VariableNames', emptyVars);

    primary = data.primaryUR;
    blocks = {
        data.lmeCoefDelta, "Delta";
        data.lmeCoefPower, "Power"};
    for bi = 1:size(blocks, 1)
        T = blocks{bi, 1};
        mclass = blocks{bi, 2};
        if isempty(T), continue; end
        T = script10_prefer_additive_lme_(T);
        for b = 1:numel(primary)
            bn = primary(b);
            if mclass == "Power" && ~(bn == "UR_1_3" || bn == "UR_3_6")
                continue;
            end
            sub = T(string(T.BandName) == bn & string(T.Term) == "Photoperiod_h", :);
            if isempty(sub), continue; end
            est = double(sub.Estimate(1));
            se = double(sub.SE(1));
            if ~isfinite(est) || ~isfinite(se), continue; end
            sig = false;
            if ismember('Significant_BH', sub.Properties.VariableNames)
                sig = logical(sub.Significant_BH(1));
            end
            rows = [rows; {mclass, bn, est, se, est - 1.96 * se, est + 1.96 * se, sig}]; %#ok<AGROW>
        end
    end
end

function T = script10_prefer_additive_lme_(T)
    if isempty(T) || ~ismember('SexModel', T.Properties.VariableNames)
        return;
    end
    sm = string(T.SexModel);
    prefer = sm == "AdditiveOnly" | contains(lower(sm), "additive");
    % Drop explicit sex-interaction model rows when AdditiveOnly available
    isSexInt = contains(lower(sm), "interaction") | contains(lower(sm), "sexint");
    if any(prefer)
        T = T(prefer, :);
    elseif any(~isSexInt)
        T = T(~isSexInt, :);
    end
end

function [widePath, tallPath] = script10_composite_fig03_(standPaths, outDirs, theme)
    imgs = cell(3, 1);
    for i = 1:3
        p = standPaths{i}; if iscell(p), p = p{1}; end
        imgs{i} = script10_trim_image_whitespace_(imread(p));
    end
    figW = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1600 900]);
    tl = tiledlayout(figW, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(tl, 1); imshow(imgs{1}); axis off;
    nexttile(tl, 2); imshow(imgs{2}); axis off;
    nexttile(tl, 3, [1 2]); imshow(imgs{3}); axis off;
    widePath = fullfile(outDirs.compositeWide, ['Fig03_Coexpression_LME' theme.ext]);
    script10_exportgraphics_(figW, widePath, theme.dpi);
    close(figW);

    figT = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1000 1400]);
    tl2 = tiledlayout(figT, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:3
        nexttile(tl2, i); imshow(imgs{i}); axis off;
    end
    tallPath = fullfile(outDirs.compositeTall, ['Fig03_Coexpression_LME' theme.ext]);
    script10_exportgraphics_(figT, tallPath, theme.dpi);
    close(figT);
end

%% Fig04 / Fig06 — Z-scored activity L12–L22 (explicit ClusterID; primary_cluster_ fallback)
function [manifest, widePath, tallPath] = script10_build_activity_figure_(data, outDirs, theme, cfg, manifest, figId, bandName, compositeStem, clusterId, writeComposite, standSubdir)
%SCRIPT10_BUILD_ACTIVITY_FIGURE_ Activity ZT grid for one cluster.
%   clusterId: ClusterID string; empty → script10_primary_cluster_
%   writeComposite: true → Wide/Tall composites (main Figs); false → Standalone only
%   standSubdir: '' → Standalone/<stem>; 'Supp' → Standalone/Supp/<stem>
    if nargin < 9, clusterId = ""; end
    if nargin < 10 || isempty(writeComposite), writeComposite = true; end
    if nargin < 11 || isempty(standSubdir), standSubdir = ''; end
    pal = theme.palette;
    if strlength(string(standSubdir)) > 0
        standDir = fullfile(outDirs.standalone, char(standSubdir), compositeStem);
    else
        standDir = fullfile(outDirs.standalone, compositeStem);
    end
    extended_period_gate_ensure_dir(standDir);
    bandName = string(bandName);
    if strlength(string(clusterId)) == 0
        clusterId = script10_primary_cluster_(data.clusterSummary, bandName);
    else
        clusterId = string(clusterId);
    end
    primaryFace = script10_cluster_face_label_(clusterId, bandName, data.clusterSummary);
    clusterCol = extended_cluster_colour(pal, 'ClusterID', clusterId, 'ClusterSummary', data.clusterSummary);
    facets = pal.coherenceFacets;
    if isfield(cfg, 'script10') && isfield(cfg.script10, 'activityFacets') && ~isempty(cfg.script10.activityFacets)
        facets = cfg.script10.activityFacets;
    elseif isfield(cfg, 'script10') && isfield(cfg.script10, 'coherenceFacets') && ~isempty(cfg.script10.coherenceFacets)
        facets = cfg.script10.coherenceFacets;
    end
    yMaxAct = script10_activity_zt_ymax_(data.activityZT, clusterId, facets);
    nMice = script10_cluster_nmice_(data.activityZT, clusterId, facets, data.nMice);

    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1600 900]);
    for fi = 1:numel(facets)
        ax = subplot(2, 3, fi); hold(ax, 'on');
        set(ax, 'Color', 'w');
        hasData = script10_plot_zt_activity_facet_(ax, data.activityZT, clusterId, facets(fi), clusterCol, pal, theme, yMaxAct);
        % Photoperiod title only — no A–F letters; per-panel n in facet helper, union n figure-level
        title(ax, char(script10_pp_label_(facets(fi))), 'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
        if ~hasData
            text(ax, 0.5, 0.5, 'No data', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontName', theme.fontName, 'FontSize', 10, 'Interpreter', 'none');
        end
        if fi == 1 || fi == 4
            ylabel(ax, 'Activity (z-scored)', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(ax, 'ZT (h)', 'FontWeight', 'bold');
        end
        set(ax, 'XLim', [0 24], 'Color', 'w');
        script10_style_axes_(ax, theme);
    end
    sgtitle(fig, sprintf('24 h activity — %s', primaryFace), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'FontSize', 14, 'Interpreter', 'none');
    script10_add_figure_n_(fig, nMice, theme);
    standPaths = script10_export_figure_(fig, fullfile(standDir, [compositeStem '_Activity']), theme, {theme.ext, '.pdf'});
    close(fig);

    widePath = '';
    tallPath = '';
    if writeComposite
        widePath = fullfile(outDirs.compositeWide, [compositeStem theme.ext]);
        copyfile(standPaths{1}, widePath, 'f');
        tallPath = fullfile(outDirs.compositeTall, [compositeStem theme.ext]);
        copyfile(standPaths{1}, tallPath, 'f');
    end

    bandCaption = script10_band_display_(bandName, 'plain');
    noteExtra = 'Profile display; panel n + union n across photoperiods.';
    if ~writeComposite
        noteExtra = [noteExtra ' Supplemental cluster export (not main Fig).'];
    end
    manifest = script10_manifest_add_(manifest, figId, 'A', standPaths{1}, 'standalone', ...
        sprintf('Cluster 24h Z-scored activity L12–L22 (%s; %s).', bandCaption, primaryFace), ...
        'Script 7', 'Tol L12/L22', noteExtra);
    if writeComposite
        manifest = script10_manifest_add_(manifest, figId, 'Composite', widePath, '16:9', ...
            sprintf('UR activity L12–L22 for %s.', primaryFace), ...
            'Script 7', 'Tol', '');
    end
end

%% Fig05 / Fig07 — 24h ZT coherence (A main) + optional LD Pre/Post R supp (BH-FDR)
function [manifest, widePath, tallPath, suppPaths] = script10_build_coherence_figure_(data, outDirs, theme, cfg, manifest, figId, bandName, compositeStem, clusterId, writeComposite, exportTransitions, standSubdir)
%SCRIPT10_BUILD_COHERENCE_FIGURE_ Panel A = 24h ZT coherence grid (main composite when writeComposite).
%   LD Pre/Post R bars are band-level Script 5 standalones under Supp_Transitions when
%   exportTransitions=true — never glued into Wide/Tall; not duplicated per cluster.
%   clusterId: ClusterID; empty → primary_cluster_
    if nargin < 9, clusterId = ""; end
    if nargin < 10 || isempty(writeComposite), writeComposite = true; end
    if nargin < 11 || isempty(exportTransitions), exportTransitions = writeComposite; end
    if nargin < 12 || isempty(standSubdir), standSubdir = ''; end
    suppPaths = strings(0, 1);
    pal = theme.palette;
    if strlength(string(standSubdir)) > 0
        standDir = fullfile(outDirs.standalone, char(standSubdir), compositeStem);
    else
        standDir = fullfile(outDirs.standalone, compositeStem);
    end
    extended_period_gate_ensure_dir(standDir);
    pp = data.ppOrder;
    ppPlot = pp(pp <= 22);
    bandName = string(bandName);

    if strlength(string(clusterId)) == 0
        clusterId = script10_primary_cluster_(data.clusterSummary, bandName);
    else
        clusterId = string(clusterId);
    end
    primaryFace = script10_cluster_face_label_(clusterId, bandName, data.clusterSummary);
    clusterCol = extended_cluster_colour(pal, 'ClusterID', clusterId, 'ClusterSummary', data.clusterSummary);
    facets = pal.coherenceFacets;
    if isfield(cfg, 'script10') && isfield(cfg.script10, 'coherenceFacets') && ~isempty(cfg.script10.coherenceFacets)
        facets = cfg.script10.coherenceFacets;
    end
    bandCol = clusterCol;
    yMax = script10_phase24_ymax_(data.phase24, clusterId, facets, pal);
    nMice = script10_cluster_nmice_phase_(data.phase24, clusterId, facets, data.nMice);

    %% A — 24 h ZT phase coherence (main composite = this panel only)
    figA = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1600 900]);
    for fi = 1:numel(facets)
        axA = subplot(2, 3, fi); hold(axA, 'on');
        set(axA, 'Color', 'w');
        hasFacet = script10_plot_coherence_zt_(axA, data.phase24, clusterId, facets(fi), ...
            bandCol, pal, theme, yMax, bandName);
        title(axA, char(script10_pp_label_(facets(fi))), 'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
        if ~hasFacet
            text(axA, 0.5, 0.5, 'No data', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontName', theme.fontName, 'FontSize', 10, 'Interpreter', 'none');
        end
        if fi == 1 || fi == 4
            ylabel(axA, 'Phase coherence R', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(axA, 'ZT (h)', 'FontWeight', 'bold');
        end
        set(axA, 'YLim', [0 yMax], 'XLim', [0 24], 'Color', 'w');
        script10_style_axes_(axA, theme);
    end
    if exportTransitions || writeComposite
        sgTxt = sprintf('A  24 h phase coherence — %s', primaryFace);
    else
        sgTxt = sprintf('24 h phase coherence — %s', primaryFace);
    end
    sgtitle(figA, sgTxt, ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'FontSize', 14, 'Interpreter', 'none');
    script10_add_figure_n_(figA, nMice, theme);
    standA = script10_export_figure_(figA, fullfile(standDir, [compositeStem '_A_CoherenceZT']), theme, {theme.ext, '.pdf'});
    close(figA);

    bandCaption = script10_band_display_(bandName, 'plain');
    noteA = 'Descriptive daily profile; panel n + union n across photoperiods.';
    if ~writeComposite
        noteA = [noteA ' Supplemental cluster export (not main Fig).'];
    end
    manifest = script10_manifest_add_(manifest, figId, 'A', standA{1}, 'standalone', ...
        sprintf('24h ZT phase coherence L12–L22 (%s; %s).', bandCaption, primaryFace), ...
        'Script 7', 'Tol cluster mean', noteA);

    widePath = '';
    tallPath = '';
    if writeComposite
        % A-only Wide/Tall (same pattern as activity figures — do not glue transition supp)
        widePath = fullfile(outDirs.compositeWide, [compositeStem theme.ext]);
        copyfile(standA{1}, widePath, 'f');
        tallPath = fullfile(outDirs.compositeTall, [compositeStem theme.ext]);
        copyfile(standA{1}, tallPath, 'f');
        manifest = script10_manifest_add_(manifest, figId, 'Composite', widePath, '16:9', ...
            sprintf('24h ZT coherence (A only) for %s. LD Pre/Post R bars are supplemental.', primaryFace), ...
            'Script 7', 'Tol', 'Main composite = panel A; LD Pre/Post R under Supp_Transitions.');
    end

    %% LD Pre/Post R — band-level (once per band; Supp_Transitions)
    if exportTransitions
        % Prefer locked coherence facets (L12–L22) when available
        if ~isempty(facets)
            ppPlot = double(facets(:))';
            ppPlot = ppPlot(ppPlot <= 22);
        end
        [manifest, suppPaths] = script10_export_transition_supp_(data, outDirs, theme, manifest, figId, bandName, ppPlot, pal);
    end
end

function [manifest, suppPaths] = script10_export_transition_supp_(data, outDirs, theme, manifest, figId, bandName, ppPlot, pal)
%SCRIPT10_EXPORT_TRANSITION_SUPP_ LD-only Pre/Post phase-coherence R bars under Supp_Transitions.
%   Exploratory display of mean R_pre / R_post; confirmatory stars from Script 5 BH-FDR (DeltaR>0).
    suppPaths = strings(0, 1);
    bandName = string(bandName);
    bandCaption = script10_band_display_(bandName, 'plain');
    figTag = char(string(figId));  % Fig05 / Fig07
    suppDir = fullfile(outDirs.standalone, 'Supp_Transitions');
    extended_period_gate_ensure_dir(suppDir);

    stem = sprintf('Supp_%s_LD_PrePost_R', figTag);
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 720 480]);
    ax = axes(fig); hold(ax, 'on');
    script10_plot_ld_prepost_r_bars_(ax, data, bandName, ppPlot, pal, theme);
    xlabel(ax, 'Photoperiod', 'FontWeight', 'bold');
    ylabel(ax, 'Phase coherence R', 'FontWeight', 'bold');
    title(ax, sprintf('LD pre/post R — %s', bandCaption), 'FontWeight', 'bold', 'Interpreter', 'none');
    script10_style_axes_(ax, theme);
    text(ax, 0.98, 0.02, '* BH-FDR \DeltaR>0 (LD)', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontName', theme.fontName, 'FontSize', 8, 'Color', [0.35 0.35 0.35], ...
        'Interpreter', 'tex');
    pathsOut = script10_export_figure_(fig, fullfile(suppDir, stem), theme, {theme.ext, '.pdf'});
    close(fig);
    suppPaths(end + 1, 1) = string(fullfile(suppDir, stem)); %#ok<AGROW>
    manifest = script10_manifest_add_(manifest, figId, 'LD_PrePost_R', pathsOut{1}, 'standalone', ...
        sprintf('LD-only Pre/Post phase-coherence R vs photoperiod (%s); BH-FDR stars (DeltaR>0).', bandCaption), ...
        'Script 5', 'Tol Pre/Post', ...
        'Confirmatory transition resync — supplemental; cite Resync_PrimaryStats_BH_FDR in text. No ridge-power / DeltaR gradient figures.');
end

function script10_plot_ld_prepost_r_bars_(ax, data, bandName, ppPlot, pal, theme)
%SCRIPT10_PLOT_LD_PREPOST_R_BARS_ Grouped Pre/Post R bars (LD) + SEM + BH-FDR brackets.
    hold(ax, 'on');
    set(ax, 'Color', 'w');
    nPP = numel(ppPlot);
    if nPP < 1
        text(ax, 0.5, 0.5, 'No photoperiods', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        return;
    end

    [meanPre, meanPost, semPre, semPost, ~] = script10_ld_prepost_r_means_(data, bandName, ppPlot);
    Y = [meanPre(:), meanPost(:)];
    if ~any(isfinite(Y(:)))
        text(ax, 0.5, 0.5, 'No LD Pre/Post R data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        set(ax, 'XTick', 1:nPP, 'XTickLabel', arrayfun(@(v) char(script10_pp_label_(v)), ppPlot, 'UniformOutput', false));
        return;
    end

    colPre = pal.base(7, :);   % Tol grey — Pre
    colPost = pal.ld;          % Tol red — Post (LD role)
    b = bar(ax, 1:nPP, Y, 0.85, 'grouped', 'EdgeColor', 'k', 'LineWidth', 0.6);
    if numel(b) >= 1
        b(1).FaceColor = colPre;
        b(1).DisplayName = 'Pre';
    end
    if numel(b) >= 2
        b(2).FaceColor = colPost;
        b(2).DisplayName = 'Post';
    end

    % Error bars at bar centres
    for bi = 1:min(2, numel(b))
        xEnds = b(bi).XEndPoints;
        if bi == 1
            y = meanPre; e = semPre;
        else
            y = meanPost; e = semPost;
        end
        for pi = 1:nPP
            if ~isfinite(y(pi)) || ~isfinite(e(pi)) || e(pi) <= 0, continue; end
            errorbar(ax, xEnds(pi), y(pi), e(pi), 'k', 'LineStyle', 'none', ...
                'LineWidth', 1.0, 'CapSize', 6, 'HandleVisibility', 'off');
        end
    end

    % BH-FDR stars / brackets for LD DeltaR>0
    fdr = data.resyncPrimaryFdr;
    sp = semPre(:); sp(~isfinite(sp)) = 0;
    sq = semPost(:); sq(~isfinite(sq)) = 0;
    yMaxData = max([meanPre(:) + sp; meanPost(:) + sq], [], 'omitnan');
    if ~isfinite(yMaxData) || yMaxData <= 0, yMaxData = 0.5; end
    yLimTop = max(yMaxData * 1.25, 0.2);
    pad = 0.04 * yLimTop;
    if ~isempty(fdr) && ismember('Significant_BH', fdr.Properties.VariableNames) && numel(b) >= 2
        xPre = b(1).XEndPoints;
        xPost = b(2).XEndPoints;
        for pi = 1:nPP
            pp = ppPlot(pi);
            idxF = fdr.Photoperiod_h == pp & string(fdr.BandName) == string(bandName) & ...
                string(fdr.TransitionType) == "LD";
            if ~any(idxF), continue; end
            sig = any(logical(fdr.Significant_BH(idxF)));
            if ~sig, continue; end
            yTop = max([meanPre(pi) + sp(pi), meanPost(pi) + sq(pi)], [], 'omitnan');
            if ~isfinite(yTop), continue; end
            y1 = yTop + pad;
            y2 = y1 + pad;
            plot(ax, [xPre(pi), xPre(pi), xPost(pi), xPost(pi)], [y1, y2, y2, y1], ...
                'k-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
            text(ax, mean([xPre(pi), xPost(pi)]), y2 + 0.5 * pad, '*', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontName', theme.fontName, 'FontSize', 14, 'FontWeight', 'bold', ...
                'HandleVisibility', 'off');
            yLimTop = max(yLimTop, y2 + 2 * pad);
        end
    end

    set(ax, 'XTick', 1:nPP, 'XTickLabel', arrayfun(@(v) char(script10_pp_label_(v)), ppPlot, 'UniformOutput', false), ...
        'XLim', [0.4, nPP + 0.6], 'YLim', [0 yLimTop]);
    lg = legend(ax, 'Location', 'northwest', 'Box', 'off');
    if ~isempty(lg)
        set(lg, 'FontName', theme.fontName, 'FontSize', 10);
    end
end

function [meanPre, meanPost, semPre, semPost, nCand] = script10_ld_prepost_r_means_(data, bandName, ppPlot)
%SCRIPT10_LD_PREPOST_R_MEANS_ Mean ± SEM R_Pre/R_Post for LD (candidate-level preferred).
%   Primary: CandidateDeltaR (PassN_Pre & PassN_Post). Fallback: DeltaR_Summary group means.
    nPP = numel(ppPlot);
    meanPre = nan(nPP, 1);
    meanPost = nan(nPP, 1);
    semPre = nan(nPP, 1);
    semPost = nan(nPP, 1);
    nCand = zeros(nPP, 1);

    CD = table();
    if isfield(data, 'candidateDeltaR') && ~isempty(data.candidateDeltaR)
        CD = data.candidateDeltaR;
    end
    useCand = ~isempty(CD) && all(ismember({'R_Pre','R_Post','Photoperiod_h','BandName','TransitionType'}, ...
        CD.Properties.VariableNames));

    if useCand
        idx = string(CD.TransitionType) == "LD" & string(CD.BandName) == string(bandName);
        if ismember('PassN_Pre', CD.Properties.VariableNames) && ismember('PassN_Post', CD.Properties.VariableNames)
            idx = idx & logical(CD.PassN_Pre) & logical(CD.PassN_Post);
        end
        C = CD(idx, :);
        for pi = 1:nPP
            rows = C(C.Photoperiod_h == ppPlot(pi), :);
            if isempty(rows), continue; end
            rp = double(rows.R_Pre); rp = rp(isfinite(rp));
            rq = double(rows.R_Post); rq = rq(isfinite(rq));
            nCand(pi) = max(numel(rp), numel(rq));
            if ~isempty(rp)
                meanPre(pi) = mean(rp, 'omitnan');
                if numel(rp) > 1
                    semPre(pi) = std(rp, 0, 'omitnan') / sqrt(numel(rp));
                else
                    semPre(pi) = 0;
                end
            end
            if ~isempty(rq)
                meanPost(pi) = mean(rq, 'omitnan');
                if numel(rq) > 1
                    semPost(pi) = std(rq, 0, 'omitnan') / sqrt(numel(rq));
                else
                    semPost(pi) = 0;
                end
            end
        end
        if any(isfinite(meanPre) | isfinite(meanPost))
            return;
        end
    end

    % Fallback: group-level DeltaR_Summary (no candidate SEM)
    DS = table();
    if isfield(data, 'deltaRSummary') && ~isempty(data.deltaRSummary)
        DS = data.deltaRSummary;
    end
    if isempty(DS) || ~all(ismember({'R_Pre','R_Post','Photoperiod_h','BandName','TransitionType'}, ...
            DS.Properties.VariableNames))
        return;
    end
    for pi = 1:nPP
        idx = DS.Photoperiod_h == ppPlot(pi) & string(DS.BandName) == string(bandName) & ...
            string(DS.TransitionType) == "LD";
        if ~any(idx), continue; end
        meanPre(pi) = double(DS.R_Pre(find(idx, 1)));
        meanPost(pi) = double(DS.R_Post(find(idx, 1)));
        semPre(pi) = 0;
        semPost(pi) = 0;
        nCand(pi) = 1;
    end
end

%% Package_Manuscript (Figs 3–7 main + Figures/Supp + tables; exclude Fig1 & Fig2)
function pkgRoot = script10_build_package_(paths, outDirs, packageFigs, packageSupp, theme, cfg, LOG) %#ok<INUSD>
    if nargin < 4 || isempty(packageSupp), packageSupp = strings(0, 1); end
    if nargin < 6 || isempty(cfg), cfg = extended_defaults(); end
    pkgRoot = fullfile(outDirs.root, 'Package_Manuscript');
    figDir = fullfile(pkgRoot, 'Figures');
    tabDir = fullfile(pkgRoot, 'Tables');
    manDir = fullfile(pkgRoot, 'Manifest');
    suppDir = fullfile(figDir, 'Supp');
    extended_period_gate_ensure_dir(figDir);
    extended_period_gate_ensure_dir(tabDir);
    extended_period_gate_ensure_dir(manDir);
    extended_period_gate_ensure_dir(suppDir);

    % Figures 3–7 composites (main only)
    for i = 1:numel(packageFigs)
        src = char(packageFigs(i));
        if isfile(src)
            copyfile(src, fullfile(figDir, basename_(src)), 'f');
        end
        [~, stem, ext] = fileparts(src);
        tallSrc = fullfile(outDirs.compositeTall, [stem ext]);
        if isfile(tallSrc)
            tallDestDir = fullfile(figDir, 'Tall_4x5');
            extended_period_gate_ensure_dir(tallDestDir);
            copyfile(tallSrc, fullfile(tallDestDir, basename_(tallSrc)), 'f');
        end
    end

    % Main standalones for Figs 3–7 (exclude Supp* folders)
    standSrcRoot = outDirs.standalone;
    standDest = fullfile(figDir, 'Standalone');
    if isfolder(standSrcRoot)
        d = dir(standSrcRoot);
        for i = 1:numel(d)
            if ~d(i).isdir || startsWith(d(i).name, '.'), continue; end
            nm = d(i).name;
            if startsWith(nm, 'Supp'), continue; end
            if startsWith(nm, 'Fig01') || startsWith(nm, 'Fig02'), continue; end
            if ~(startsWith(nm, 'Fig03') || startsWith(nm, 'Fig04') || startsWith(nm, 'Fig05') || ...
                    startsWith(nm, 'Fig06') || startsWith(nm, 'Fig07'))
                continue;
            end
            dest = fullfile(standDest, nm);
            extended_period_gate_ensure_dir(fileparts(dest));
            copyfile(fullfile(d(i).folder, nm), dest, 'f');
        end
    end

    % Supplemental figures → Figures/Supp/ (LD Pre/Post R + C02 grids)
    for i = 1:numel(packageSupp)
        src = char(packageSupp(i));
        if isfolder(src)
            [~, leaf] = fileparts(src);
            dest = fullfile(suppDir, leaf);
            extended_period_gate_ensure_dir(fileparts(dest));
            if isfolder(dest)
                try, rmdir(dest, 's'); catch, end %#ok<CTCH>
            end
            copyfile(src, dest, 'f');
        elseif isfile(src)
            copyfile(src, fullfile(suppDir, basename_(src)), 'f');
        else
            % Stem folder may hold PNG+PDF siblings — copy matching files from parent
            parent = fileparts(src);
            stem = basename_(src);
            if isfolder(parent)
                matches = dir(fullfile(parent, [stem '.*']));
                for mi = 1:numel(matches)
                    if matches(mi).isdir, continue; end
                    copyfile(fullfile(matches(mi).folder, matches(mi).name), ...
                        fullfile(suppDir, matches(mi).name), 'f');
                end
                % Also copy whole Supp_Transitions / Supp subtrees once
            end
        end
    end
    % Ensure full Supp_Transitions and Supp trees are in package
    for subName = {'Supp_Transitions', 'Supp'}
        srcTree = fullfile(outDirs.standalone, subName{1});
        if isfolder(srcTree)
            destTree = fullfile(suppDir, subName{1});
            extended_period_gate_ensure_dir(fileparts(destTree));
            if isfolder(destTree)
                try, rmdir(destTree, 's'); catch, end %#ok<CTCH>
            end
            copyfile(srcTree, destTree, 'f');
        end
    end

    % Supplemental tables (Scripts 4–7 + optional 9/11) → Package + Script10/Tables/Supplemental
    destDirs = struct();
    destDirs.packageTables = tabDir;
    destDirs.script10Tables = fullfile(outDirs.root, 'Tables', 'Supplemental');
    extended_script10_export_supplemental_tables(paths, destDirs, cfg, LOG);

    % Legacy flat CSV aliases at Tables/ root (backward-compatible handoff paths)
    script10_copy_priority_csv_aliases_(tabDir);

    % Manifest copies
    if isfile(outDirs.manifest)
        copyfile(outDirs.manifest, fullfile(manDir, 'Manifest.xlsx'), 'f');
    end
    mdSrc = fullfile(outDirs.root, 'FIGURE_MANIFEST.md');
    if isfile(mdSrc)
        copyfile(mdSrc, fullfile(manDir, 'FIGURE_MANIFEST.md'), 'f');
    end
    csvSrc = fullfile(outDirs.root, 'FigureManifest.csv');
    if isfile(csvSrc)
        copyfile(csvSrc, fullfile(manDir, 'FigureManifest.csv'), 'f');
    end

    script10_write_package_readme_(pkgRoot, theme);
end

function script10_write_package_readme_(pkgRoot, theme) %#ok<INUSD>
    fid = fopen(fullfile(pkgRoot, 'README.md'), 'w', 'n', 'UTF-8');
    if fid < 0, return; end
    fprintf(fid, '# Package_Manuscript - Script 10 handoff\n\n');
    fprintf(fid, 'Zip this folder for collaborator / manuscript handoff.\n\n');
    fprintf(fid, '## Included\n\n');
    fprintf(fid, '- **Figures/** - Fig03-Fig07 main composites (Wide + optional Tall/Standalone)\n');
    fprintf(fid, '- **Figures/Supp/** - supplemental panels (not renumbered mains):\n');
    fprintf(fid, '  - `Supp_Transitions/` - Fig05/07 **LD Pre/Post R** bars (Script 5 BH-FDR stars)\n');
    fprintf(fid, '  - `Supp/` - UR_3_6 **C02** activity + coherence twins\n');
    fprintf(fid, '- **Tables/** - full supplemental export (`Workbooks/`, `CSV/`, `README_SUPPLEMENTAL_TABLES.md`)\n');
    fprintf(fid, '  - Methods/QC (Script 4 CarryForward), co-expression + LME (Script 6),\n');
    fprintf(fid, '    transition resync + LL projected (Script 5), profiles (Script 7),\n');
    fprintf(fid, '    Script 9 manifest when present (Script 11 tables stay in Script11_DominantPeriod_*)\n');
    fprintf(fid, '  - Flat CSV aliases at `Tables/*.csv` retained for legacy paths\n');
    fprintf(fid, '- **Manifest/** - figure order, captions, source scripts\n\n');
    fprintf(fid, '## Excluded (still in full Script10 tree)\n\n');
    fprintf(fid, '- **Fig01** - circadian characteristics placeholder (collaborator external)\n');
    fprintf(fid, '- **Fig02** - RAW/HSub scalograms + CarryForward retention QC\n\n');
    fprintf(fid, '## Scientific boundaries (do not conflate)\n\n');
    fprintf(fid, '| Figure | Inference source | Claim family |\n');
    fprintf(fid, '|--------|------------------|--------------|\n');
    fprintf(fid, '| **Fig03** A-B | Descriptive means | Exploratory co-expression gradients |\n');
    fprintf(fid, '| **Fig03C** | Script 6 LME Photoperiod_h beta+/-CI (BH-FDR) | **Confirmatory co-expression** |\n');
    fprintf(fid, '| **Figs 04 & 06** | Script 7 profiles (C01) | Primary-cluster 24h activity display |\n');
    fprintf(fid, '| **Figs 05 & 07** | Script 7 24h ZT coherence (**A only** on main) | Descriptive daily phase organisation |\n');
    fprintf(fid, '| **Supp LD Pre/Post R (Fig05/07)** | Script 5 BH-FDR DeltaR>0 (LD bars) | **Confirmatory transition resync** (cite tables; panels supplemental) |\n');
    fprintf(fid, '| **Supp UR36 C02** | Script 7 profiles (ClusterRank 2) | Extra cluster face; not main Fig 8/9 |\n\n');
    fprintf(fid, 'Do not conflate Fig 3 LME with Script 5 transition stats.\n');
    fprintf(fid, 'Main Figs 5/7 composites are coherence grids only; LD Pre/Post R live under Figures/Supp/Supp_Transitions/.\n');
    fprintf(fid, 'Sex-stratified panels remain Script 9 (supplemental).\n');
    fclose(fid);

    fid2 = fopen(fullfile(pkgRoot, 'Manifest', 'PACKAGE_README.md'), 'w');
    if fid2 > 0
        fprintf(fid2, 'See `../README.md` for included vs excluded content and claim boundaries.\n');
        fclose(fid2);
    end
end

function script10_copy_if_exists_(src, dest)
    if isfile(src)
        extended_period_gate_ensure_dir(fileparts(dest));
        copyfile(src, dest, 'f');
    end
end

function script10_copy_priority_csv_aliases_(tabDir)
%SCRIPT10_COPY_PRIORITY_CSV_ALIASES Flat aliases at Tables/ root for legacy handoff paths.
    csvSub = fullfile(tabDir, 'CSV');
    if ~isfolder(csvSub), return; end
    aliases = { ...
        'CR_UR_Pairs_Summary.csv'; ...
        'AbsolutePower_Summary.csv'; ...
        'LME_Coef_Delta_BH_FDR.csv'; ...
        'LME_Coef_Power_BH_FDR.csv'; ...
        'BinnedCoherence.csv'; ...
        'RidgePowerStats_BH_FDR.csv'; ...
        'Resync_PrimaryStats_BH_FDR.csv'; ...
        'TransitionEffect_vs_Photoperiod_Summary.csv'; ...
        'ClusterSummary.csv'; ...
        'ActivityComponent_24h.csv'; ...
        'PhaseCoherence_24h.csv'; ...
        'ClusterMembership.csv'};
    for i = 1:numel(aliases)
        src = fullfile(csvSub, aliases{i});
        if isfile(src)
            copyfile(src, fullfile(tabDir, aliases{i}), 'f');
        end
    end
end

%% Shared plotting helpers (adapted from Script 8; local — not exported)
function cid = script10_primary_cluster_(CS, bandName)
%SCRIPT10_PRIMARY_CLUSTER_ ClusterRank 1 for band (legacy C01 path).
    cid = script10_cluster_by_rank_(CS, bandName, 1);
end

function cid = script10_cluster_by_rank_(CS, bandName, rank)
%SCRIPT10_CLUSTER_BY_RANK_ Resolve ClusterID from ClusterSummary by ClusterRank.
    cid = "";
    if isempty(CS) || ~ismember('BandName', CS.Properties.VariableNames), return; end
    sub = CS(string(CS.BandName) == string(bandName), :);
    if isempty(sub), return; end
    rank = double(rank);
    if ismember('ClusterRank', sub.Properties.VariableNames)
        hit = sub(double(sub.ClusterRank) == rank, :);
        if ~isempty(hit)
            cid = string(hit.ClusterID(1));
            return;
        end
        % Fallback: sort by rank and take ordinal if exact rank missing
        sub = sortrows(sub, 'ClusterRank', 'ascend');
        if rank >= 1 && rank <= height(sub)
            cid = string(sub.ClusterID(rank));
        end
        return;
    end
    if ismember('CandidateCount', sub.Properties.VariableNames)
        sub = sortrows(sub, 'CandidateCount', 'descend');
        if rank >= 1 && rank <= height(sub)
            cid = string(sub.ClusterID(rank));
        end
    end
end

function cid = script10_resolve_manuscript_cluster_(CS, bandName, rank, cfg, LOG)
%SCRIPT10_RESOLVE_MANUSCRIPT_CLUSTER_ ClusterRank lock + optional period-window warn.
    if nargin < 5, LOG = -1; end
    cid = script10_cluster_by_rank_(CS, bandName, rank);
    if strlength(string(cid)) == 0
        script10_log_(LOG, 'No cluster for %s rank %g', char(string(bandName)), double(rank));
        return;
    end
    if nargin < 4 || isempty(cfg) || ~isfield(cfg, 'script10') || ...
            ~isfield(cfg.script10, 'manuscriptClusters')
        return;
    end
    mc = cfg.script10.manuscriptClusters;
    bandKey = matlab.lang.makeValidName(char(string(bandName)));
    if ~isfield(mc, bandKey), return; end
    spec = mc.(bandKey);
    if ~isfield(spec, 'expectPeriodRanks') || ~isfield(spec, 'expectPeriod_h') || isempty(spec.expectPeriod_h)
        return;
    end
    er = double(spec.expectPeriodRanks(:));
    ep = double(spec.expectPeriod_h);
    idx = find(er == double(rank), 1);
    if isempty(idx) || idx > size(ep, 1), return; end
    loExp = ep(idx, 1); hiExp = ep(idx, 2);
    if ~all(ismember({'ClusterID', 'PeriodLow_h', 'PeriodHigh_h'}, CS.Properties.VariableNames)), return; end
    row = CS(string(CS.ClusterID) == string(cid), :);
    if isempty(row), return; end
    lo = double(row.PeriodLow_h(1)); hi = double(row.PeriodHigh_h(1));
    if ~(isfinite(lo) && isfinite(hi)), return; end
    tol = 0.15;  % hours — soft check against Script 7 window drift
    if abs(lo - loExp) > tol || abs(hi - hiExp) > tol
        msg = sprintf(['Cluster period window drift for %s rank %g (%s): got %.2f–%.2f h, ' ...
            'expected ~%.2f–%.2f h (cfg.script10.manuscriptClusters).'], ...
            char(string(bandName)), double(rank), char(string(cid)), lo, hi, loExp, hiExp);
        warning('extended_script10:ClusterPeriodCheck', '%s', msg);
        script10_log_(LOG, '%s', msg);
    end
end

function ranks = script10_manuscript_supp_ranks_(cfg, bandName)
    ranks = [];
    if ~isfield(cfg, 'script10') || ~isfield(cfg.script10, 'manuscriptClusters'), return; end
    mc = cfg.script10.manuscriptClusters;
    bandKey = matlab.lang.makeValidName(char(string(bandName)));
    if ~isfield(mc, bandKey), return; end
    spec = mc.(bandKey);
    if isfield(spec, 'suppRanks') && ~isempty(spec.suppRanks)
        ranks = double(spec.suppRanks(:))';
    end
end

function rank = script10_manuscript_main_rank_(cfg, bandName)
%SCRIPT10_MANUSCRIPT_MAIN_RANK_ First mainRanks entry (default 1 = C01).
    rank = 1;
    if ~isfield(cfg, 'script10') || ~isfield(cfg.script10, 'manuscriptClusters'), return; end
    mc = cfg.script10.manuscriptClusters;
    bandKey = matlab.lang.makeValidName(char(string(bandName)));
    if ~isfield(mc, bandKey), return; end
    spec = mc.(bandKey);
    if isfield(spec, 'mainRanks') && ~isempty(spec.mainRanks)
        rank = double(spec.mainRanks(1));
    end
end

function lbl = script10_cluster_short_label_(clusterID, CS)
    lbl = char(string(clusterID));
    if isempty(CS) || strlength(string(clusterID)) == 0, return; end
    if ~ismember('ClusterID', CS.Properties.VariableNames), return; end
    sub = CS(string(CS.ClusterID) == string(clusterID), :);
    if isempty(sub), return; end
    if ismember('ClusterRank', sub.Properties.VariableNames)
        lbl = sprintf('C%02d', double(sub.ClusterRank(1)));
    end
end

function lbl = script10_cluster_face_label_(clusterID, bandName, CS)
    band = script10_band_display_(bandName, 'plain');
    short = script10_cluster_short_label_(clusterID, CS);
    if strlength(string(short)) == 0
        lbl = band; return;
    end
    rangeTxt = '';
    if ~isempty(CS) && ismember('ClusterID', CS.Properties.VariableNames)
        sub = CS(string(CS.ClusterID) == string(clusterID), :);
        if ~isempty(sub) && all(ismember({'PeriodLow_h', 'PeriodHigh_h'}, sub.Properties.VariableNames))
            lo = double(sub.PeriodLow_h(1)); hi = double(sub.PeriodHigh_h(1));
            if isfinite(lo) && isfinite(hi)
                rangeTxt = sprintf(', %.1f–%.1f h', lo, hi);
            end
        end
    end
    lbl = sprintf('%s (%s%s)', band, short, rangeTxt);
end

function yLim = script10_activity_zt_ymax_(Act, clusterID, facets)
    yLim = [-0.5 1.5];
    if isempty(Act), return; end
    if strlength(string(clusterID)) == 0, return; end
    needed = {'ClusterID', 'Photoperiod_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames)), return; end
    A = Act(string(Act.ClusterID) == string(clusterID) & ismember(Act.Photoperiod_h, facets), :);
    if isempty(A), return; end
    y = double(A.Activity_zscored); y = y(isfinite(y));
    if isempty(y), return; end
    pad = 0.12 * max(range(y), 0.5);
    yLim = [min(y) - pad, max(y) + pad];
end

function hasData = script10_plot_zt_activity_facet_(ax, Act, clusterID, photo, meanCol, pal, theme, yLim)
    if nargin < 8, yLim = []; end
    hasData = false;
    if isempty(Act) || strlength(string(clusterID)) == 0, return; end
    needed = {'ClusterID', 'Photoperiod_h', 'SignalID', 'ZTBinCenter_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames)), return; end
    A = Act(string(Act.ClusterID) == string(clusterID) & Act.Photoperiod_h == photo, :);
    if isempty(A), return; end
    hold(ax, 'on'); set(ax, 'Color', 'w');
    if ~ismember('File', A.Properties.VariableNames)
        A.File = repmat("", height(A), 1);
    end
    mouseKey = string(A.File) + "|" + string(A.SignalID);
    sigs = unique(mouseKey, 'stable');
    mouseZT = {}; mouseY = {};
    for s = 1:numel(sigs)
        As = A(mouseKey == sigs(s), :);
        As = sortrows(As, 'ZTBinCenter_h');
        zt = double(As.ZTBinCenter_h); y = double(As.Activity_zscored);
        keep = zt >= 0 & zt <= 24 & isfinite(y);
        if ~any(keep), continue; end
        plot(ax, zt(keep), y(keep), '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 0.85, 'HandleVisibility', 'off');
        mouseZT{end+1,1} = zt(keep); %#ok<AGROW>
        mouseY{end+1,1} = y(keep); %#ok<AGROW>
    end
    if ~isempty(mouseZT)
        hasData = true;
        edges = 0:0.5:24; centers = edges(1:end-1) + diff(edges)/2;
        meanY = nan(size(centers));
        for i = 1:numel(centers)
            vals = nan(numel(mouseZT), 1);
            for m = 1:numel(mouseZT)
                idx = mouseZT{m} >= edges(i) & mouseZT{m} < edges(i+1);
                if i == numel(centers)
                    idx = mouseZT{m} >= edges(i) & mouseZT{m} <= edges(i+1);
                end
                if any(idx), vals(m) = mean(mouseY{m}(idx), 'omitnan'); end
            end
            meanY(i) = mean(vals, 'omitnan');
        end
        plot(ax, centers, meanY, '-', 'Color', meanCol, 'LineWidth', 2.6);
        script10_add_panel_n_(ax, numel(mouseZT), theme);
    end
    if ~isempty(yLim)
        ylim(ax, yLim); script10_shade_zt_ld_(ax, double(photo), pal, yLim);
    else
        yl = ylim(ax);
        if diff(yl) < 1e-6, yl = [-0.5 1.5]; ylim(ax, yl); end
        script10_shade_zt_ld_(ax, double(photo), pal, yl);
    end
    xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    photoH = double(photo);
    if photoH > 0 && photoH < 24
        xline(ax, photoH, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
    xlim(ax, [0 24]);
    script10_style_axes_(ax, theme);
end

function script10_shade_zt_ld_(ax, photoH, pal, yl) %#ok<INUSD>
%SCRIPT10_SHADE_ZT_LD_ Dark-phase shade + lights-off line (collaborator / Script 9 style).
    if nargin < 4 || isempty(yl), yl = ylim(ax); end
    if numel(yl) < 2 || ~all(isfinite(yl)) || yl(2) <= yl(1)
        yl = [-1 1];
    end
    photoH = double(photoH);
    if ~isfinite(photoH) || photoH >= 24 || photoH <= 0, return; end
    h = patch(ax, [photoH 24 24 photoH], [yl(1) yl(1) yl(2) yl(2)], [0.88 0.88 0.88], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    if ~isempty(h), uistack(h, 'bottom'); end
    xline(ax, photoH, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

function yMax = script10_coherence_ymax_(Bin, facets, bandName, pal)
    yMax = pal.coherenceYMax;
    if isempty(Bin) || ~ismember('R', Bin.Properties.VariableNames), return; end
    B = Bin(ismember(Bin.Photoperiod_h, facets) & string(Bin.BandName) == string(bandName), :);
    if isempty(B), return; end
    mx = max(double(B.R), [], 'omitnan');
    if isfinite(mx) && mx > 0
        yMax = max(ceil(mx * 1.05 * 20) / 20, 0.2);
    end
end

function hasData = script10_plot_coherence_facet_(ax, Bin, photoperiod_h, bandName, pal, theme)
    hasData = false;
    if isempty(Bin), return; end
    B = Bin(Bin.Photoperiod_h == photoperiod_h & string(Bin.BandName) == bandName, :);
    if isempty(B), return; end
    hold(ax, 'on'); set(ax, 'Color', 'w');
    for tr = ["DL", "LD"]
        Bt = B(string(B.TransitionType) == tr, :);
        if isempty(Bt), continue; end
        hasData = true;
        Bt = sortrows(Bt, 'RelBinCenter_h');
        col = script10_transition_colour_(pal, tr);
        x = Bt.RelBinCenter_h; y = Bt.R;
        plot(ax, x, y, '-o', 'Color', col, 'LineWidth', 2.0, 'MarkerSize', 4, 'MarkerFaceColor', col);
        if ismember('N_PhaseObs', Bt.Properties.VariableNames)
            se = 1.96 ./ sqrt(max(Bt.N_PhaseObs, 1));
            fill(ax, [x; flipud(x)], [y - se; flipud(y + se)], col, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        end
        if x(end) > 3
            script10_direct_line_label_(ax, x(end), y(end), char(tr), col, theme, false);
        end
    end
    xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.3, 'HandleVisibility', 'off');
    script10_style_axes_(ax, theme);
end

function script10_plot_empty_coherence_note_(ax, photoperiod_h, theme)
    axis(ax, [0 1 0 1]); axis(ax, 'off');
    text(ax, 0.5, 0.55, sprintf('No data for %s', char(script10_pp_label_(photoperiod_h))), ...
        'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 10);
end

%% Theme / IO / manifest
function theme = script10_theme_(cfg)
    pal = extended_tol_bright_palette();
    dpi = double(cfg.plot.saveDpi);
    if ~isscalar(dpi) || ~isfinite(dpi) || dpi <= 0, dpi = 150; end
    ext = char(string(cfg.plot.figExt));
    if strlength(string(ext)) == 0
        ext = '.png';
    elseif ~startsWith(ext, '.'), ext = ['.' ext];
    end
    theme = struct('palette', pal, 'fontName', pal.fontName, ...
        'dpi', round(dpi), 'ext', ext, 'plotMode', string(cfg.plotMode));
end

function [outDirs, modeLabel] = script10_make_output_dirs_(cohortRoot, plotMode)
    pm = lower(strtrim(char(string(plotMode))));
    if strcmp(pm, 'publication')
        modeLabel = "Publication";
    else
        modeLabel = "Development";
        if ~strcmp(pm, 'development') && strlength(string(pm)) > 0
            warning('extended_script10:UnknownPlotMode', ...
                'Unrecognized plotMode "%s"; using Development folder.', pm);
        end
    end
    outRoot = fullfile(cohortRoot, ['Script10_ManuscriptFigures_' char(modeLabel)]);
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

function manifest = script10_manifest_init_()
    manifest = table('Size', [0 8], ...
        'VariableTypes', repmat({'string'}, 1, 8), ...
        'VariableNames', {'Figure', 'Panel', 'File', 'Aspect', 'Caption', 'SourceScript', 'ColourKey', 'Notes'});
end

function manifest = script10_manifest_add_(manifest, figId, panelId, filePath, aspect, caption, sourceScript, colourKey, notes)
    if nargin < 9, notes = ""; end
    manifest = [manifest; {string(figId), string(panelId), string(filePath), string(aspect), ...
        string(caption), string(sourceScript), string(colourKey), string(notes)}]; %#ok<AGROW>
end

function script10_write_manifest_(manifest, xlsxPath)
    writetable(manifest, xlsxPath, 'Sheet', 'Manifest');
    writetable(script10_colour_key_table_(), xlsxPath, 'Sheet', 'ColourKey');
    % Captions sheet for narrative lock
    caps = manifest(manifest.Panel == "Composite" | manifest.Panel == "Placeholder" | manifest.Panel == "A", ...
        {'Figure','Panel','Caption','SourceScript','Notes'});
    writetable(caps, xlsxPath, 'Sheet', 'Captions');
end

function script10_write_figure_manifest_md_(outRoot, manifest)
    fid = fopen(fullfile(outRoot, 'FIGURE_MANIFEST.md'), 'w', 'n', 'UTF-8');
    if fid < 0, return; end
    fprintf(fid, '# Script 10 - Locked manuscript figure set\n\n');
    fprintf(fid, '| Fig | Content | Source | Confirmatory vs exploratory |\n');
    fprintf(fid, '|-----|---------|--------|-----------------------------|\n');
    fprintf(fid, '| 1 | Circadian characteristics - **placeholder** (collaborator) | External | Slot only |\n');
    fprintf(fid, '| 2 | Multiscale RAW + HSub residual Female\\|Male + CarryForward E | Scripts 3/4 | Display / QC |\n');
    fprintf(fid, '| 3 | A abs power; B descriptive Delta; **C LME forest** | Script 6 | A-B exploratory; **C confirmatory co-expression** |\n');
    fprintf(fid, '| 4 | UR_1_3 Z-scored activity L12-L22 (**C01**) | Script 7 | Profile display |\n');
    fprintf(fid, '| 5 | UR_1_3 **24h ZT coherence (A only on main)** | Script 7 | Descriptive daily phase |\n');
    fprintf(fid, '| 6 | UR_3_6 activity twin of Fig 4 (**C01**) | Script 7 | Profile display |\n');
    fprintf(fid, '| 7 | UR_3_6 coherence twin of Fig 5 (**A only on main**) | Script 7 | Descriptive daily phase |\n\n');
    fprintf(fid, '### Supplemental (not main storyboard / not renumbered mains)\n\n');
    fprintf(fid, '| Asset | Content | Source |\n');
    fprintf(fid, '|-------|---------|--------|\n');
    fprintf(fid, '| `Standalone/Supp_Transitions/Supp_Fig05_LD_PrePost_R` | UR_1_3 LD Pre/Post R bars + BH-FDR | Script 5 |\n');
    fprintf(fid, '| `Standalone/Supp_Transitions/Supp_Fig07_LD_PrePost_R` | UR_3_6 LD Pre/Post R bars + BH-FDR | Script 5 |\n');
    fprintf(fid, '| `Standalone/Supp/Supp_Activity_UR36_C02` | UR_3_6 ClusterRank 2 activity twin | Script 7 |\n');
    fprintf(fid, '| `Standalone/Supp/Supp_Coherence_UR36_C02` | UR_3_6 ClusterRank 2 coherence twin | Script 7 |\n\n');
    fprintf(fid, 'Cluster lock: `cfg.script10.manuscriptClusters` (ClusterRank; optional period-window warn).\n');
    fprintf(fid, 'Main Figs 4-7 = C01 only; LD Pre/Post R and C02 are Package `Figures/Supp/`.\n');
    fprintf(fid, 'Cite Script 5 `Resync_PrimaryStats_BH_FDR` in text for transition claims (panels supplemental).\n');
    fprintf(fid, 'Ridge-power / DeltaR gradient figures are not exported (tables remain in Package).\n');
    fprintf(fid, 'Sex after Fig 2 -> Script 9 supplemental.\n\n');
    fprintf(fid, '## Manifest rows\n\n');
    for i = 1:height(manifest)
        fprintf(fid, '- **%s / %s** - %s (%s)\n', ...
            char(manifest.Figure(i)), char(manifest.Panel(i)), ...
            char(manifest.Caption(i)), char(manifest.SourceScript(i)));
    end
    fclose(fid);
end

function script10_write_figure_manifest_csv_(outRoot, manifest)
    try
        writetable(manifest, fullfile(outRoot, 'FigureManifest.csv'));
    catch
    end
end

function T = script10_colour_key_table_()
    roles = {
        'UR 1–3 C01 mean', '#66CCEE', 'Activity + coherence twins; Script 11 pooled violin';
        'UR 3–6 C01 mean', '#AA3377', 'Activity + coherence twins; Script 11 pooled violin';
        'UR 3–6 C02 mean', '#EE6677', 'Supp activity + coherence twins; Script 11 pooled violin (cluster mean, not LD transition bar)';
        'Individual mice', '#B3B3B3', 'Light grey ZT traces (all cluster panels)';
        'LD Pre R', '#BBBBBB', 'Supp LD Pre/Post R bars';
        'LD Post R', '#EE6677', 'Supp LD Pre/Post R bars (Tol LD)';
        'DL (lights-on)', '#4477AA', 'Transition traces only (not cluster mean on ZT grids)';
        'LD (lights-off)', '#EE6677', 'Transition colour role';
        'Female', '#228833', 'Sex-stratified panels only';
        'Male', '#CCBB44', 'Sex-stratified panels only';
        'UR 1–3 h band (legacy)', '#66CCEE', 'Band-level panels (Fig 3, LME forest)';
        'UR 3–6 h band (legacy)', '#AA3377', 'Band-level panels (Fig 3, LME forest)';
        'CR 20–28 h', '#BBBBBB', 'CR traces';
        'BH-FDR star', 'black *', 'Script 5 LD DeltaR>0 (Supp) or Script 6 (Fig 3C)';
        'Scalograms', 'Jet', 'Fig02 RAW + HSub residual'};
    T = cell2table(roles, 'VariableNames', {'Role', 'Colour', 'Usage'});
end

function script10_save_storyboard_(compositeFiles, pdfPath)
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

function script10_export_tall_copy_(tallPath, widePath, tallDir)
    if isfile(tallPath), return; end
    if isfile(widePath)
        copyfile(widePath, fullfile(tallDir, ['TALL_' basename_(widePath)]), 'f');
    end
end

function b = basename_(p)
    [~, b, e] = fileparts(p);
    b = [b e];
end

function img = script10_trim_image_whitespace_(img)
    if isempty(img), return; end
    if ndims(img) == 3
        mask = any(img < 250, 3);
    else
        mask = img < 250;
    end
    [rows, cols] = find(mask);
    if isempty(rows) || isempty(cols), return; end
    pad = 6;
    r1 = max(min(rows) - pad, 1); r2 = min(max(rows) + pad, size(img, 1));
    c1 = max(min(cols) - pad, 1); c2 = min(max(cols) + pad, size(img, 2));
    img = img(r1:r2, c1:c2, :);
end

function img = script10_crop_hsub_residual_(imgPath)
    img = imread(imgPath);
    img = img(:, round(size(img, 2) / 2) + 1:end, :);
end

function script10_exportgraphics_(fig, outFile, dpi)
    if nargin < 3 || isempty(dpi), dpi = 150; end
    dpi = round(double(dpi(1)));
    if ~isfinite(dpi) || dpi <= 0, dpi = 150; end
    try
        [outDir, baseName, ext] = fileparts(outFile);
        if strlength(string(ext)) == 0
            ext = '.png'; outFile = fullfile(outDir, [baseName ext]);
        end
        extLower = lower(ext);
        basePath = fullfile(outDir, baseName);
        if ismember(extLower, {'.pdf', '.eps', '.emf'})
            exportgraphics(fig, outFile, 'ContentType', 'vector'); return;
        end
        if script10_try_exportgraphics_raster_(fig, outFile, dpi), return; end
        resFlag = sprintf('-r%d', dpi);
        produced = '';
        switch extLower
            case {'.jpg', '.jpeg'}
                print(fig, basePath, '-djpeg', resFlag, '-noui'); produced = [basePath '.jpg'];
            case '.png'
                print(fig, basePath, '-dpng', resFlag, '-noui'); produced = [basePath '.png'];
            otherwise
                script10_saveas_fallback_(fig, outFile, basePath, extLower); return;
        end
        if isfile(outFile), return; end
        if strlength(string(produced)) > 0 && isfile(produced) && ~strcmpi(produced, outFile)
            copyfile(produced, outFile, 'f'); return;
        end
        script10_saveas_fallback_(fig, outFile, basePath, extLower);
    catch ME
        warning('Script10:ExportFailed', 'Could not export figure to %s: %s', outFile, ME.message);
    end
end

function script10_saveas_fallback_(fig, outFile, basePath, extLower)
    try
        if ismember(extLower, {'.jpg', '.jpeg'})
            saveas(fig, basePath, 'jpg');
            produced = [basePath '.jpg'];
            if ~strcmpi(produced, outFile) && isfile(produced), copyfile(produced, outFile, 'f'); end
        elseif ismember(extLower, {'.png'})
            saveas(fig, basePath, 'png');
            produced = [basePath '.png'];
            if ~strcmpi(produced, outFile) && isfile(produced), copyfile(produced, outFile, 'f'); end
        else
            saveas(fig, outFile);
        end
    catch ME
        warning('Script10:ExportFailed', 'Could not export figure to %s: %s', outFile, ME.message);
    end
end

function ok = script10_try_exportgraphics_raster_(fig, outFile, dpi)
    ok = false;
    for k = 1:4
        switch k
            case 1, resVal = dpi;
            case 2, resVal = uint16(dpi);
            case 3, resVal = int32(dpi);
            otherwise, resVal = [];
        end
        try
            if isempty(resVal), exportgraphics(fig, outFile);
            else, exportgraphics(fig, outFile, 'Resolution', resVal);
            end
            ok = true; return;
        catch
        end
    end
end

function outPaths = script10_export_figure_(fig, basePath, theme, formats)
    outPaths = cell(numel(formats), 1);
    extended_period_gate_ensure_dir(fileparts(basePath));
    [bpDir, bpName, bpExt] = fileparts(basePath);
    if strlength(string(bpExt)) > 0, basePath = fullfile(bpDir, bpName); end
    dpi = double(theme.dpi);
    if ~isscalar(dpi) || ~isfinite(dpi) || dpi <= 0, dpi = 150; end
    dpi = round(dpi);
    for i = 1:numel(formats)
        ext = char(formats{i});
        if ~startsWith(ext, '.'), ext = ['.' ext]; end
        outFile = [basePath ext];
        script10_exportgraphics_(fig, outFile, dpi);
        outPaths{i} = outFile;
    end
end

function T = script10_read_sheet_(path, sheet)
    try
        T = readtable(path, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    catch
        T = table();
    end
end

function ppOrder = script10_pp_order_(T)
    if isempty(T) || ~ismember('Photoperiod_h', T.Properties.VariableNames)
        ppOrder = [12 14 16 18 20 22 24]; return;
    end
    ppOrder = unique(double(T.Photoperiod_h));
    ppOrder = sort(ppOrder(isfinite(ppOrder)));
end

function n = script10_total_mice_(sexBalance)
    n = 24;
    if isempty(sexBalance) || ~ismember('N_Mice', sexBalance.Properties.VariableNames), return; end
    if ismember('Sex', sexBalance.Properties.VariableNames)
        G = findgroups(sexBalance.Sex);
        mx = splitapply(@max, sexBalance.N_Mice, G);
        n = sum(mx);
    else
        n = max(sexBalance.N_Mice);
    end
end

function lbl = script10_pp_label_(pp)
    lbl = "L" + string(round(double(pp)));
end

function lbl = script10_band_display_(bandKey, mode)
    if nargin < 2, mode = 'tex'; end
    key = char(string(bandKey));
    switch key
        case 'UR_1_3', texLbl = 'UR_{1–3}'; plainLbl = 'UR 1–3 h';
        case 'UR_3_6', texLbl = 'UR_{3–6}'; plainLbl = 'UR 3–6 h';
        case 'UR_6_9', texLbl = 'UR_{6–9}'; plainLbl = 'UR 6–9 h';
        case 'UR_9_12', texLbl = 'UR_{9–12}'; plainLbl = 'UR 9–12 h';
        case 'UR_12_18', texLbl = 'UR_{12–18}'; plainLbl = 'UR 12–18 h';
        case 'CR_20_28', texLbl = 'CR_{20–28}'; plainLbl = 'CR 20–28 h';
        otherwise, texLbl = strrep(key, '_', '\_'); plainLbl = strrep(key, '_', ' ');
    end
    if strcmpi(mode, 'plain'), lbl = plainLbl; else, lbl = texLbl; end
end

function rgb = script10_band_colour_(pal, bandName)
    key = char(string(bandName));
    if isKey(pal.band, key), rgb = pal.band(key); else, rgb = pal.base(7, :); end
end

function rgb = script10_transition_colour_(pal, transitionType)
    t = upper(string(transitionType));
    if t == "DL", rgb = pal.dl; elseif t == "LD", rgb = pal.ld; else, rgb = pal.base(7, :); end
end

function script10_style_axes_(ax, theme)
    set(ax, 'Box', 'off', 'TickDir', theme.palette.tickDir, 'FontName', theme.fontName, ...
        'LineWidth', theme.palette.axesLineWidth, 'FontSize', 11);
    if isprop(ax, 'XLabel') && ~isempty(ax.XLabel.String)
        ax.XLabel.FontWeight = 'bold'; ax.XLabel.FontSize = 12;
    end
    if isprop(ax, 'YLabel') && ~isempty(ax.YLabel.String)
        ax.YLabel.FontWeight = 'bold'; ax.YLabel.FontSize = 12;
    end
end

function script10_panel_label_(ax, labelChar, theme)
    text(ax, 0.02, 0.98, labelChar, 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 16, 'FontName', theme.fontName, 'VerticalAlignment', 'top');
end

function script10_direct_line_label_(ax, x, y, txt, colour, theme, useTex)
    if nargin < 7, useTex = false; end
    interp = 'none'; if useTex, interp = 'tex'; end
    text(ax, x, y, txt, 'Color', colour, 'FontName', theme.fontName, 'FontWeight', 'bold', ...
        'FontSize', 10, 'Interpreter', interp);
end

function script10_add_n_annotation_(ax, nMice, theme)
    text(ax, 0.98, 0.02, sprintf('n = %d mice', nMice), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontName', theme.fontName, 'FontSize', 9);
end

function script10_add_panel_n_(ax, nMice, theme)
%SCRIPT10_ADD_PANEL_N_ Per-facet mouse count (cluster members plotted in this panel).
    if ~(isfinite(nMice) && nMice > 0), return; end
    text(ax, 0.98, 0.98, sprintf('n = %d', round(nMice)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontName', theme.fontName, 'FontSize', 9, 'Color', [0.25 0.25 0.25], ...
        'Interpreter', 'none', 'HandleVisibility', 'off');
end

function script10_add_figure_n_(fig, nMice, theme)
%SCRIPT10_ADD_FIGURE_N_ Bottom-right union n across all facets (not per panel).
    if ~(isfinite(nMice) && nMice > 0), return; end
    annotation(fig, 'textbox', [0.55 0.005 0.44 0.035], ...
        'String', sprintf('n = %d mice (across photoperiods)', round(nMice)), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontName', theme.fontName, 'FontSize', 11, 'FontWeight', 'bold', ...
        'Color', [0.15 0.15 0.15], 'Interpreter', 'none');
end

function n = script10_cluster_nmice_(Act, clusterID, facets, fallbackN)
%SCRIPT10_CLUSTER_NMICE_ Unique mice union across primary-cluster activity facets.
%   Identity is animal-stable across photoperiods: prefer MouseID, else SignalID.
%   Do NOT key on File|SignalID — File differs per photoperiod and inflates n.
%   Per-panel counts can be lower (cluster membership only when a validated ridge is present).
    n = fallbackN;
    if isempty(Act) || strlength(string(clusterID)) == 0, return; end
    if ~all(ismember({'ClusterID', 'Photoperiod_h'}, Act.Properties.VariableNames)), return; end
    A = Act(string(Act.ClusterID) == string(clusterID) & ismember(Act.Photoperiod_h, facets), :);
    if isempty(A), return; end
    keys = script10_unique_mouse_keys_(A);
    if ~isempty(keys)
        n = numel(keys);
    end
end

function n = script10_cluster_nmice_phase_(Phase24, clusterID, facets, fallbackN)
%SCRIPT10_CLUSTER_NMICE_PHASE_ Unique mice for coherence figure n annotation.
%   Prefer animal ID (MouseID / SignalID); else max(N_Mice) on aggregated bins.
    n = fallbackN;
    if isempty(Phase24) || strlength(string(clusterID)) == 0, return; end
    B = Phase24(ismember(Phase24.Photoperiod_h, facets), :);
    if ismember('ClusterID', B.Properties.VariableNames)
        B = B(string(B.ClusterID) == string(clusterID), :);
    end
    if isempty(B), return; end
    keys = script10_unique_mouse_keys_(B);
    if ~isempty(keys)
        n = numel(keys);
        return;
    end
    if ismember('N_Mice', B.Properties.VariableNames)
        nAgg = max(double(B.N_Mice), [], 'omitnan');
        if isfinite(nAgg) && nAgg > 0
            n = nAgg;
        end
    end
end

function keys = script10_unique_mouse_keys_(T)
%SCRIPT10_UNIQUE_MOUSE_KEYS_ Stable animal IDs (no File / photoperiod in key).
    keys = strings(0, 1);
    if isempty(T), return; end
    if ismember('MouseID', T.Properties.VariableNames)
        keys = unique(string(T.MouseID));
        keys = keys(strlength(keys) > 0);
        if ~isempty(keys), return; end
    end
    if ismember('SignalID', T.Properties.VariableNames)
        keys = unique(string(T.SignalID));
        keys = keys(strlength(keys) > 0);
    end
end

function yMax = script10_phase24_ymax_(Phase24, clusterID, photos, pal)
    yMax = 0.8;
    if isfield(pal, 'coherenceYMax') && isfinite(pal.coherenceYMax)
        yMax = max(0.2, double(pal.coherenceYMax));
    end
    if isempty(Phase24) || ~ismember('R', Phase24.Properties.VariableNames), return; end
    B = Phase24(ismember(Phase24.Photoperiod_h, photos), :);
    if ismember('ClusterID', B.Properties.VariableNames) && strlength(string(clusterID)) > 0
        B = B(string(B.ClusterID) == string(clusterID), :);
    end
    if isempty(B), return; end
    mx = max(double(B.R), [], 'omitnan');
    if isfinite(mx) && mx > 0
        yMax = max(ceil(mx * 1.08 * 20) / 20, 0.2);
    end
end

function has = script10_plot_coherence_zt_(ax, Phase24, clusterID, photo, lineCol, pal, theme, yMax, bandName) %#ok<INUSL>
    has = false;
    if isempty(Phase24), return; end
    needed = {'Photoperiod_h', 'ZTBinCenter_h', 'R'};
    if ~all(ismember(needed, Phase24.Properties.VariableNames)), return; end
    B = Phase24(Phase24.Photoperiod_h == photo, :);
    if ismember('ClusterID', B.Properties.VariableNames) && strlength(string(clusterID)) > 0
        B = B(string(B.ClusterID) == string(clusterID), :);
    elseif ismember('BandName', B.Properties.VariableNames)
        B = B(string(B.BandName) == string(bandName), :);
    end
    if isempty(B), return; end
    has = true;
    hold(ax, 'on'); set(ax, 'Color', 'w');
    script10_shade_zt_ld_(ax, double(photo), pal, [0 yMax]);
    B = sortrows(B, 'ZTBinCenter_h');
    if ismember('SignalID', B.Properties.VariableNames)
        sigs = unique(string(B.SignalID), 'stable');
        for s = 1:numel(sigs)
            Bs = B(string(B.SignalID) == sigs(s), :);
            plot(ax, Bs.ZTBinCenter_h, Bs.R, '-', 'Color', [0.75 0.75 0.75], ...
                'LineWidth', 0.7, 'HandleVisibility', 'off');
        end
        script10_add_panel_n_(ax, numel(sigs), theme);
    elseif ismember('N_Mice', B.Properties.VariableNames)
        nAgg = max(double(B.N_Mice), [], 'omitnan');
        if isfinite(nAgg) && nAgg > 0
            script10_add_panel_n_(ax, nAgg, theme);
        end
    end
    G = groupsummary(B, 'ZTBinCenter_h', 'mean', 'R');
    if ismember('mean_R', G.Properties.VariableNames)
        plot(ax, G.ZTBinCenter_h, G.mean_R, '-', 'Color', lineCol, 'LineWidth', 2.4);
    else
        plot(ax, B.ZTBinCenter_h, B.R, '-', 'Color', lineCol, 'LineWidth', 2.4);
    end
end

function script10_yline_zero_(ax)
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
end

function script10_log_(LOG, fmt, varargin)
    if ~isempty(LOG) && LOG > 0
        fprintf(LOG, ['[' datestr(now, 'yyyy-mm-dd HH:MM:SS') '] ' fmt '\n'], varargin{:});
    end
end

function script10_fclose_(LOG)
    if ~isempty(LOG) && LOG > 0, fclose(LOG); end
end
