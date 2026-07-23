function out = extended_script9_supplementary_figures_run(cohortRoot, cfg)
%EXTENDED_SCRIPT9_SUPPLEMENTARY_FIGURES_RUN Publication supplementary figures.
%
%   Outputs (per primary cluster + extras):
%     S01 Sex-stratified CR–UR delta
%     S02/S03 UR 1–3 C01: 24h coherence 2x3 + activity 2x3 (L12–L22)
%     S04/S05 UR 3–6 C01: 24h coherence 2x3 + activity 2x3 (L12–L22)
%     S06+   Other clusters: activity 2x3 each
%     Optional: coexpression RAW copy
%
%   Export uses short temp paths then copyfile (avoids long-path / .jpeg saveas failures).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    cfg = extended_apply_plot_cfg(cfg);
    theme = script9_theme_(cfg);

    if nargin < 1 || isempty(cohortRoot)
        cohortRoot = uigetdir(pwd, 'Select cohort results folder (e.g. C57_LP)');
        if isequal(cohortRoot, 0)
            error('extended_script9_supplementary_figures_run:NoRoot', 'No cohort folder selected.');
        end
    end

    paths = extended_script8_resolve_paths(cohortRoot);
    need = intersect(paths.missing, ["profilesXlsx", "lmeDescriptiveXlsx", "lmeInferenceXlsx"]);
    if ~isempty(need)
        error('extended_script9_supplementary_figures_run:MissingInputs', ...
            'Missing required inputs: %s', strjoin(need, ', '));
    end

    outDirs = script9_make_output_dirs_(paths.cohortRoot, cfg.plotMode);
    logPath = fullfile(outDirs.logs, sprintf('Script9_Supplementary_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() script9_fclose_(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 9: Supplementary figures ===\n');
    fprintf('Cohort:  %s\n', paths.cohortRoot);
    fprintf('Output:  %s\n', outDirs.root);
    fprintf('Mode:    %s | dpi=%g | ext=%s\n', cfg.plotMode, theme.dpi, theme.ext);
    script9_log_(LOG, 'Cohort: %s', paths.cohortRoot);

    data = script9_load_data_(paths, cfg);
    manifest = script9_manifest_init_();
    storyboard = strings(0, 1);
    nOk = 0;

    %% S01 — Sex delta
    try
        [manifest, p] = script9_build_sex_(data, outDirs, theme, manifest);
        if strlength(p) > 0
            storyboard(end+1,1) = p; %#ok<AGROW>
            nOk = nOk + 1;
        end
        script9_log_(LOG, 'S01 sex complete: %s', p);
    catch ME
        warning('Script9:Sex', 'S01 sex failed: %s', ME.message);
        script9_log_(LOG, 'S01 FAILED: %s', ME.message);
    end

    %% Primary-cluster coherence + activity (2x3 L12–L22)
    facets = theme.palette.coherenceFacets;
    nBands = min(2, numel(data.primaryUR));
    for bi = 1:nBands
        bn = data.primaryUR(bi);
        cid = script9_primary_cluster_(data.clusterSummary, bn);
        if strlength(cid) == 0
            warning('Script9:NoCluster', 'No primary cluster for %s', char(bn));
            continue;
        end
        face = script9_cluster_face_label_(cid, bn, data.clusterSummary);
        stem = script9_band_file_stem_(bn);
        cidShort = script9_cluster_short_label_(cid, data.clusterSummary);
        fileStem = sprintf('%s_%s', stem, cidShort);

        try
            fig = script9_make_coherence_fig_(data, bn, cid, facets, theme, face);
            outP = script9_export_jpeg_(fig, fullfile(outDirs.figures, ['Supp_ClusterCoherence_' fileStem]), theme);
            close(fig);
            if strlength(outP) == 0
                error('Export returned empty path');
            end
            storyboard(end+1,1) = outP; %#ok<AGROW>
            nOk = nOk + 1;
            manifest = script9_manifest_add_(manifest, 'Supp', ['Coh_' fileStem], outP, ...
                sprintf('24h phase coherence L12–L22 — %s', face));
            script9_log_(LOG, 'Coherence %s: %s', fileStem, outP);
        catch ME
            warning('Script9:Coherence', 'Coherence %s failed: %s', fileStem, ME.message);
            script9_log_(LOG, 'Coherence %s FAILED: %s', fileStem, ME.message);
        end

        try
            fig = script9_make_activity_fig_(data, cid, facets, theme, face);
            outP = script9_export_jpeg_(fig, fullfile(outDirs.figures, ['Supp_ClusterActivity_' fileStem]), theme);
            close(fig);
            if strlength(outP) == 0
                error('Export returned empty path');
            end
            storyboard(end+1,1) = outP; %#ok<AGROW>
            nOk = nOk + 1;
            manifest = script9_manifest_add_(manifest, 'Supp', ['Act_' fileStem], outP, ...
                sprintf('24h ZT activity L12–L22 — %s', face));
            script9_log_(LOG, 'Activity %s: %s', fileStem, outP);
        catch ME
            warning('Script9:Activity', 'Activity %s failed: %s', fileStem, ME.message);
            script9_log_(LOG, 'Activity %s FAILED: %s', fileStem, ME.message);
        end
    end

    %% Other clusters — activity only
    CS = data.clusterSummary;
    primaryIds = strings(0, 1);
    for bi = 1:nBands
        primaryIds(end+1,1) = script9_primary_cluster_(CS, data.primaryUR(bi)); %#ok<AGROW>
    end
    if ~isempty(CS) && ismember('ClusterID', CS.Properties.VariableNames)
        for ci = 1:height(CS)
            cid = string(CS.ClusterID(ci));
            if strlength(cid) == 0 || any(cid == primaryIds)
                continue;
            end
            bn = "";
            if ismember('BandName', CS.Properties.VariableNames)
                bn = string(CS.BandName(ci));
            end
            face = script9_cluster_face_label_(cid, bn, CS);
            cidSafe = regexprep(char(cid), '[^A-Za-z0-9_]+', '_');
            try
                fig = script9_make_activity_fig_(data, cid, facets, theme, face);
                outP = script9_export_jpeg_(fig, fullfile(outDirs.figures, ['Supp_Activity_' cidSafe]), theme);
                close(fig);
                if strlength(outP) == 0
                    error('Export returned empty path');
                end
                storyboard(end+1,1) = outP; %#ok<AGROW>
                nOk = nOk + 1;
                manifest = script9_manifest_add_(manifest, 'Supp', ['ActOther_' cidSafe], outP, ...
                    sprintf('24h activity L12–L22 — non-primary %s', face));
                script9_log_(LOG, 'Other activity %s: %s', cidSafe, outP);
            catch ME
                warning('Script9:OtherAct', 'Other activity %s failed: %s', cidSafe, ME.message);
                script9_log_(LOG, 'Other activity %s FAILED: %s', cidSafe, ME.message);
            end
        end
    end

    %% Optional coexpression copy
    coexpPath = fullfile(paths.script3Root, 'Figures', 'Coexpression_CR_UR_RAW.jpeg');
    if isfile(coexpPath)
        dest = fullfile(outDirs.figures, ['Supp_Coexpression_RAW' theme.ext]);
        try
            copyfile(coexpPath, dest, 'f');
            storyboard(end+1,1) = string(dest); %#ok<AGROW>
            nOk = nOk + 1;
            manifest = script9_manifest_add_(manifest, 'Supp', 'Coexpression', dest, ...
                'Core exploratory RAW co-expression (Script 3).');
        catch
        end
    end

    script9_write_manifest_(manifest, outDirs.manifest);
    script9_save_storyboard_(storyboard, outDirs.storyboard);

    out = struct();
    out.cohortRoot = paths.cohortRoot;
    out.outRoot = outDirs.root;
    out.manifest = outDirs.manifest;
    out.storyboard = outDirs.storyboard;
    out.nFigures = nOk;
    out.logPath = logPath;

    fprintf('\nExtended Script 9 complete (%d figure(s)).\n', nOk);
    fprintf('  Figures:    %s\n', outDirs.figures);
    fprintf('  Storyboard: %s\n', outDirs.storyboard);
    fprintf('  Manifest:   %s\n', outDirs.manifest);
    if nOk == 0
        warning('Script9:NoFigures', 'No supplementary figures were written. Check log: %s', logPath);
    end
end

%% ------------------------------------------------------------------------
function data = script9_load_data_(paths, cfg)
    data = struct();
    data.pairSummary = script9_read_sheet_(paths.lmeDescriptiveXlsx, 'CR_UR_Pairs_Summary');
    data.pairBySex = script9_read_sheet_(paths.lmeDescriptiveXlsx, 'CR_UR_Pairs_Summary_BySex');
    data.phase24 = script9_read_sheet_(paths.profilesXlsx, 'PhaseCoherence_24h');
    data.clusterSummary = script9_read_sheet_(paths.profilesXlsx, 'ClusterSummary');
    data.activityZT = script9_read_sheet_(paths.profilesXlsx, 'ActivityComponent_24h');
    data.primaryUR = string(cfg.bands.primaryUR);
    data.ppOrder = script9_pp_order_(data.pairSummary);
end

function [manifest, outPath] = script9_build_sex_(data, outDirs, theme, manifest)
    pal = theme.palette;
    pp = data.ppOrder;
    ppLabels = arrayfun(@(x) char(script9_pp_label_(x)), pp, 'UniformOutput', false);
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 460]);
    for bi = 1:numel(data.primaryUR)
        ax = subplot(1, 2, bi); hold(ax, 'on');
        bn = data.primaryUR(bi);
        if ~isempty(data.pairBySex)
            pooled = data.pairSummary(string(data.pairSummary.UR_Band) == bn & string(data.pairSummary.Phase) == "All", :);
            pooled = sortrows(pooled, 'Photoperiod_h');
            if ~isempty(pooled)
                plot(ax, 1:height(pooled), pooled.Mean_Delta_log10, '--', 'Color', pal.pooled, 'LineWidth', 1.2);
            end
            for sx = ["Female", "Male"]
                sub = data.pairBySex(string(data.pairBySex.UR_Band) == bn & string(data.pairBySex.Phase) == "All" & string(data.pairBySex.Sex) == sx, :);
                sub = sortrows(sub, 'Photoperiod_h');
                if isempty(sub), continue; end
                if sx == "Female", col = pal.female; else, col = pal.male; end
                plot(ax, 1:height(sub), sub.Mean_Delta_log10, '-o', 'Color', col, 'LineWidth', 2.2, 'MarkerFaceColor', col);
                text(ax, height(sub) + 0.05, sub.Mean_Delta_log10(end), char(sx), ...
                    'Color', col, 'FontName', theme.fontName, 'FontSize', 10, 'Interpreter', 'none');
            end
        end
        yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
        set(ax, 'XTick', 1:numel(pp), 'XTickLabel', ppLabels);
        xlabel(ax, 'Photoperiod', 'FontWeight', 'bold');
        if bi == 1
            ylabel(ax, '\Delta(UR - CR) log_{10}', 'FontWeight', 'bold');
        end
        title(ax, ['Sex-stratified \Delta — ' script9_band_display_(bn, 'tex')], ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
        script9_style_axes_(ax, theme);
    end
    outPath = script9_export_jpeg_(fig, fullfile(outDirs.figures, 'Supp_Fig03D_Sex'), theme);
    close(fig);
    if strlength(outPath) > 0
        manifest = script9_manifest_add_(manifest, 'Supp', 'Sex_delta', outPath, ...
            'Sex-stratified CR–UR \Delta (off main Fig03).');
    end
end

function fig = script9_make_coherence_fig_(data, bandName, clusterID, facets, theme, faceLabel)
    pal = theme.palette;
    yMax = script9_coherence_ymax_(data.phase24, clusterID, facets, pal);
    bandCol = script9_band_colour_(pal, bandName);
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1600 900]);
    for fi = 1:numel(facets)
        ax = subplot(2, 3, fi); hold(ax, 'on'); set(ax, 'Color', 'w');
        has = script9_plot_coherence_zt_(ax, data.phase24, clusterID, facets(fi), bandCol, pal, theme, yMax, bandName);
        title(ax, char(script9_pp_label_(facets(fi))), 'FontWeight', 'bold', 'Interpreter', 'none');
        if ~has
            text(ax, 0.5, 0.5, sprintf('No coherence for %s', char(script9_pp_label_(facets(fi)))), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontName', theme.fontName, 'Interpreter', 'none');
        end
        if fi == 1 || fi == 4
            ylabel(ax, 'Phase coherence R', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(ax, 'ZT (h)', 'FontWeight', 'bold');
        end
        set(ax, 'XLim', [0 24], 'YLim', [0 yMax], 'Color', 'w');
        script9_style_axes_(ax, theme);
        script9_panel_label_(ax, char('A' + fi - 1), theme);
    end
    sgtitle(fig, sprintf('24h phase coherence — %s', faceLabel), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
end

function fig = script9_make_activity_fig_(data, clusterID, facets, theme, faceLabel)
    pal = theme.palette;
    yLim = script9_activity_ymax_(data.activityZT, clusterID, facets);
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [60 60 1600 900]);
    for fi = 1:numel(facets)
        ax = subplot(2, 3, fi); hold(ax, 'on'); set(ax, 'Color', 'w');
        has = script9_plot_activity_zt_(ax, data.activityZT, clusterID, facets(fi), pal, theme, yLim);
        title(ax, char(script9_pp_label_(facets(fi))), 'FontWeight', 'bold', 'Interpreter', 'none');
        if ~has
            text(ax, 0.5, 0.5, sprintf('No activity for %s', char(script9_pp_label_(facets(fi)))), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontName', theme.fontName, 'Interpreter', 'none');
        end
        if fi == 1 || fi == 4
            ylabel(ax, 'Activity (z-scored)', 'FontWeight', 'bold');
        end
        if fi >= 4
            xlabel(ax, 'ZT (h)', 'FontWeight', 'bold');
        end
        set(ax, 'XLim', [0 24], 'YLim', yLim, 'Color', 'w');
        script9_style_axes_(ax, theme);
        script9_panel_label_(ax, char('A' + fi - 1), theme);
    end
    sgtitle(fig, sprintf('24h activity — %s', faceLabel), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
end

function has = script9_plot_coherence_zt_(ax, Phase24, clusterID, photo, lineCol, pal, theme, yMax, bandName) %#ok<INUSL>
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
    script9_shade_ld_(ax, photo, [0 yMax]);
    B = sortrows(B, 'ZTBinCenter_h');
    if ismember('SignalID', B.Properties.VariableNames)
        sigs = unique(string(B.SignalID), 'stable');
        for s = 1:numel(sigs)
            Bs = B(string(B.SignalID) == sigs(s), :);
            plot(ax, Bs.ZTBinCenter_h, Bs.R, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7, 'HandleVisibility', 'off');
        end
    end
    G = groupsummary(B, 'ZTBinCenter_h', 'mean', 'R');
    if ismember('mean_R', G.Properties.VariableNames)
        plot(ax, G.ZTBinCenter_h, G.mean_R, '-o', 'Color', lineCol, 'LineWidth', 2.4, ...
            'MarkerSize', 4, 'MarkerFaceColor', lineCol);
    end
    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
end

function has = script9_plot_activity_zt_(ax, Act, clusterID, photo, pal, theme, yLim) %#ok<INUSL>
    has = false;
    if isempty(Act) || strlength(string(clusterID)) == 0, return; end
    needed = {'ClusterID', 'Photoperiod_h', 'SignalID', 'ZTBinCenter_h', 'Activity_zscored'};
    if ~all(ismember(needed, Act.Properties.VariableNames)), return; end
    A = Act(string(Act.ClusterID) == string(clusterID) & Act.Photoperiod_h == photo, :);
    if isempty(A), return; end
    has = true;
    script9_shade_ld_(ax, photo, yLim);
    if ~ismember('File', A.Properties.VariableNames)
        A.File = repmat("", height(A), 1);
    end
    mouseKey = string(A.File) + "|" + string(A.SignalID);
    sigs = unique(mouseKey, 'stable');
    allZ = [];
    allY = [];
    for s = 1:numel(sigs)
        As = sortrows(A(mouseKey == sigs(s), :), 'ZTBinCenter_h');
        z = double(As.ZTBinCenter_h);
        y = double(As.Activity_zscored);
        keep = isfinite(z) & isfinite(y);
        if ~any(keep), continue; end
        plot(ax, z(keep), y(keep), '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 0.85, 'HandleVisibility', 'off');
        allZ = [allZ; z(keep)]; %#ok<AGROW>
        allY = [allY; y(keep)]; %#ok<AGROW>
    end
    if ~isempty(allZ)
        edges = 0:0.5:24;
        centers = edges(1:end-1) + diff(edges)/2;
        meanY = nan(size(centers));
        for i = 1:numel(centers)
            idx = allZ >= edges(i) & allZ < edges(i+1);
            if i == numel(centers)
                idx = allZ >= edges(i) & allZ <= edges(i+1);
            end
            if any(idx)
                meanY(i) = mean(allY(idx), 'omitnan');
            end
        end
        plot(ax, centers, meanY, '-', 'Color', pal.ld, 'LineWidth', 2.6);
    end
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
end

function script9_shade_ld_(ax, photoH, yl)
    photoH = double(photoH);
    if ~isfinite(photoH) || photoH >= 24
        return;
    end
    if nargin < 3 || isempty(yl) || numel(yl) < 2 || ~all(isfinite(yl))
        yl = [-3 3];
    end
    patch(ax, [photoH 24 24 photoH], [yl(1) yl(1) yl(2) yl(2)], [0.88 0.88 0.88], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    xline(ax, photoH, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

function yMax = script9_coherence_ymax_(Phase24, clusterID, photos, pal)
    yMax = 1.0;
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
        yMax = max(ceil(mx * 1.05 * 20) / 20, 0.2);
    end
end

function yLim = script9_activity_ymax_(Act, clusterID, facets)
    yLim = [-3 3];
    if isempty(Act) || strlength(string(clusterID)) == 0, return; end
    if ~ismember('Activity_zscored', Act.Properties.VariableNames), return; end
    A = Act(string(Act.ClusterID) == string(clusterID) & ismember(Act.Photoperiod_h, facets), :);
    if isempty(A), return; end
    y = double(A.Activity_zscored);
    y = y(isfinite(y));
    if isempty(y), return; end
    m = max(abs(y));
    if isfinite(m) && m > 0
        pad = ceil(m * 1.1 * 2) / 2;
        yLim = [-pad pad];
    end
end

function outPath = script9_export_jpeg_(fig, destBaseNoExt, theme)
% Export via short temp path then copyfile (robust vs long paths / .jpeg saveas).
    outPath = '';
    extended_period_gate_ensure_dir(fileparts(destBaseNoExt));
    dpi = round(double(theme.dpi));
    if ~isscalar(dpi) || ~isfinite(dpi) || dpi <= 0, dpi = 150; end
    dpi = max(72, dpi);

    destBase = char(destBaseNoExt);
    [~, ~, e0] = fileparts(destBase);
    if strlength(string(e0)) > 0
        destBase = destBase(1:end-numel(e0));
    end
    wantExt = lower(char(string(theme.ext)));
    if ~ismember(wantExt, {'.jpg', '.jpeg', '.png'})
        wantExt = '.jpeg';
    end
    finalPath = [destBase wantExt];

    tmpDir = fullfile(tempdir, 'BioRhythms_Script9');
    extended_period_gate_ensure_dir(tmpDir);
    tmpBase = fullfile(tmpDir, sprintf('s9_%s', datestr(now, 'HHMMSSFFF')));

    drawnow;
    ok = false;

    try
        print(fig, tmpBase, '-djpeg', sprintf('-r%d', dpi), '-noui');
        produced = [tmpBase '.jpg'];
        if isfile(produced)
            copyfile(produced, finalPath, 'f');
            ok = isfile(finalPath);
            try, delete(produced); catch, end
        end
    catch
    end

    if ~ok
        try
            tmpJpg = [tmpBase '_eg.jpg'];
            exportgraphics(fig, tmpJpg);
            if isfile(tmpJpg)
                copyfile(tmpJpg, finalPath, 'f');
                ok = isfile(finalPath);
                try, delete(tmpJpg); catch, end
            end
        catch
        end
    end

    if ~ok
        try
            set(fig, 'Visible', 'on');
            drawnow;
            fr = getframe(fig);
            set(fig, 'Visible', 'off');
            if ~isempty(fr) && isfield(fr, 'cdata') && ~isempty(fr.cdata)
                imwrite(fr.cdata, finalPath, 'jpg', 'Quality', 95);
                ok = isfile(finalPath);
            end
        catch
            try, set(fig, 'Visible', 'off'); catch, end
        end
    end

    if ok
        outPath = finalPath;
    else
        warning('Script9:ExportFailed', 'Could not export figure to %s', finalPath);
    end
end

function cid = script9_primary_cluster_(CS, bandName)
    cid = "";
    if isempty(CS) || ~ismember('BandName', CS.Properties.VariableNames), return; end
    sub = CS(string(CS.BandName) == string(bandName), :);
    if isempty(sub), return; end
    if ismember('ClusterRank', sub.Properties.VariableNames)
        sub = sortrows(sub, 'ClusterRank', 'ascend');
    elseif ismember('CandidateCount', sub.Properties.VariableNames)
        sub = sortrows(sub, 'CandidateCount', 'descend');
    end
    cid = string(sub.ClusterID(1));
end

function lbl = script9_cluster_short_label_(clusterID, CS)
    lbl = 'C01';
    if ~isempty(CS) && ismember('ClusterID', CS.Properties.VariableNames) && ismember('ClusterRank', CS.Properties.VariableNames)
        hit = CS(string(CS.ClusterID) == string(clusterID), :);
        if ~isempty(hit)
            lbl = sprintf('C%02d', double(hit.ClusterRank(1)));
            return;
        end
    end
    tok = regexp(char(clusterID), 'C(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        lbl = sprintf('C%02d', str2double(tok{1}));
    end
end

function lbl = script9_cluster_face_label_(clusterID, bandName, CS)
    bandDisp = script9_band_display_(bandName, 'plain');
    shortC = script9_cluster_short_label_(clusterID, CS);
    periodTxt = '';
    if ~isempty(CS) && ismember('ClusterID', CS.Properties.VariableNames)
        hit = CS(string(CS.ClusterID) == string(clusterID), :);
        if ~isempty(hit) && all(ismember({'PeriodLow_h','PeriodHigh_h'}, hit.Properties.VariableNames))
            periodTxt = sprintf(', %.1f–%.1f h', double(hit.PeriodLow_h(1)), double(hit.PeriodHigh_h(1)));
        end
    end
    lbl = sprintf('%s (%s%s)', bandDisp, shortC, periodTxt);
end

function stem = script9_band_file_stem_(bandName)
    bn = char(string(bandName));
    if contains(bn, '1_3') || contains(bn, '1-3')
        stem = 'UR13';
    elseif contains(bn, '3_6') || contains(bn, '3-6')
        stem = 'UR36';
    else
        stem = regexprep(bn, '[^A-Za-z0-9]+', '');
    end
end

function lbl = script9_pp_label_(pp)
    lbl = sprintf('L%d', round(double(pp)));
end

function lbl = script9_band_display_(bandKey, mode)
    if nargin < 2, mode = 'plain'; end
    k = char(string(bandKey));
    switch k
        case {'UR_1_3', 'UR1_3'}
            if strcmpi(mode, 'tex'), lbl = 'UR_{1–3}'; else, lbl = 'UR 1–3 h'; end
        case {'UR_3_6', 'UR3_6'}
            if strcmpi(mode, 'tex'), lbl = 'UR_{3–6}'; else, lbl = 'UR 3–6 h'; end
        case {'CR_20_28'}
            if strcmpi(mode, 'tex'), lbl = 'CR_{20–28}'; else, lbl = 'CR 20–28 h'; end
        otherwise
            lbl = k;
    end
end

function rgb = script9_band_colour_(pal, bandName)
    bn = char(string(bandName));
    if isfield(pal, 'band') && isa(pal.band, 'containers.Map') && isKey(pal.band, bn)
        rgb = pal.band(bn);
    elseif contains(bn, '1_3')
        rgb = pal.base(2, :);
    elseif contains(bn, '3_6')
        rgb = pal.base(6, :);
    else
        rgb = pal.base(1, :);
    end
end

function script9_style_axes_(ax, theme)
    set(ax, 'FontName', theme.fontName, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'off');
    ax.XLabel.FontWeight = 'bold'; ax.XLabel.FontSize = 12;
    ax.YLabel.FontWeight = 'bold'; ax.YLabel.FontSize = 12;
end

function script9_panel_label_(ax, labelChar, theme)
    text(ax, 0.02, 0.98, labelChar, 'Units', 'normalized', 'FontName', theme.fontName, ...
        'FontWeight', 'bold', 'FontSize', 16, 'VerticalAlignment', 'top', 'Interpreter', 'none');
end

function theme = script9_theme_(cfg)
    pal = extended_tol_bright_palette();
    dpi = double(cfg.plot.saveDpi);
    if ~isscalar(dpi) || ~isfinite(dpi) || dpi <= 0, dpi = 150; end
    ext = char(string(cfg.plot.figExt));
    if ~startsWith(ext, '.'), ext = ['.' ext]; end
    theme = struct('palette', pal, 'fontName', pal.fontName, 'dpi', dpi, 'ext', ext, ...
        'plotMode', string(cfg.plotMode));
end

function outDirs = script9_make_output_dirs_(cohortRoot, plotMode)
    modeLabel = "Development";
    if strcmpi(string(plotMode), "publication")
        modeLabel = "Publication";
    end
    outRoot = fullfile(cohortRoot, char("Script9_SupplementaryFigures_" + modeLabel));
    outDirs = struct( ...
        'root', outRoot, ...
        'figures', fullfile(outRoot, 'Figures'), ...
        'logs', fullfile(outRoot, 'Logs'), ...
        'manifest', fullfile(outRoot, 'Manifest.xlsx'), ...
        'storyboard', fullfile(outRoot, 'Storyboard_Supplementary.pdf'));
    extended_period_gate_ensure_dir(outDirs.figures);
    extended_period_gate_ensure_dir(outDirs.logs);
end

function manifest = script9_manifest_init_()
    manifest = table('Size', [0 5], ...
        'VariableTypes', {'string','string','string','string','string'}, ...
        'VariableNames', {'FigureID','PanelID','FilePath','Caption','Notes'});
end

function manifest = script9_manifest_add_(manifest, figId, panelId, filePath, caption)
    row = {string(figId), string(panelId), string(filePath), string(caption), ""};
    manifest = [manifest; cell2table(row, 'VariableNames', manifest.Properties.VariableNames)]; %#ok<AGROW>
end

function script9_write_manifest_(manifest, xlsxPath)
    try
        writetable(manifest, xlsxPath, 'Sheet', 'Manifest');
    catch ME
        warning('Script9:Manifest', 'Could not write manifest: %s', ME.message);
    end
end

function script9_save_storyboard_(files, pdfPath)
    files = files(strlength(files) > 0);
    files = files(isfile(files));
    if isempty(files), return; end
    for i = 1:numel(files)
        try
            img = imread(files(i));
            fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
                'Position', [40 40 max(size(img,2), 400) max(size(img,1), 300)]);
            ax = axes(fig, 'Position', [0 0 1 1]); %#ok<LAXES>
            imshow(img, 'Parent', ax); axis(ax, 'off');
            if i == 1
                exportgraphics(fig, pdfPath, 'ContentType', 'image');
            else
                exportgraphics(fig, pdfPath, 'ContentType', 'image', 'Append', true);
            end
            close(fig);
        catch
        end
    end
end

function T = script9_read_sheet_(path, sheet)
    try
        T = readtable(path, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    catch
        T = table();
    end
end

function ppOrder = script9_pp_order_(T)
    if isempty(T) || ~ismember('Photoperiod_h', T.Properties.VariableNames)
        ppOrder = [12 14 16 18 20 22 24];
        return;
    end
    ppOrder = unique(double(T.Photoperiod_h));
    ppOrder = sort(ppOrder(isfinite(ppOrder)));
end

function script9_log_(LOG, fmt, varargin)
    if LOG <= 0, return; end
    fprintf(LOG, '[%s] %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), sprintf(fmt, varargin{:}));
end

function script9_fclose_(LOG)
    if LOG > 0, try, fclose(LOG); catch, end; end
end
