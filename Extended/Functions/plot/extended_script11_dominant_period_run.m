function out = extended_script11_dominant_period_run(cohortRoot, cfg)
%EXTENDED_SCRIPT11_DOMINANT_PERIOD_RUN Dominant cluster period tables + supp violins.
%
%   Inputs (read-only; no WP_TS):
%     Script 7  SelectedValidatedUR_*Profiles/Tables/*_Output.xlsx
%               ClusterSummary + ClusterMembership (RawPeriod_h per candidate)
%
%   Per mouse: mean + median of candidate RawPeriod_h within each locked cluster.
%   Cohort: mean + median across mice (of those per-mouse summaries).
%
%   Outputs under {cohortRoot}/Script11_DominantPeriod_{Publication|Development}/
%     Tables/DominantPeriod_ClusterSummary.csv
%     Tables/DominantPeriod_PerMouse.csv
%     Tables/DominantPeriod_Output.xlsx
%     Figures/Supp_DominantPeriod_ClusterPeriod_3Cluster.{jpeg|png}

    SCRIPT_NAME = 'extended_script11_dominant_period_run';
    SCRIPT_VERSION = '1.1';

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    if ~isfield(cfg, 'plotMode') || isempty(cfg.plotMode)
        cfg.plotMode = 'development';
    end
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);
    cfg = script11_fill_cfg_(cfg);

    if nargin < 1 || isempty(cohortRoot)
        cohortRoot = uigetdir(pwd, 'Select cohort results folder (e.g. C57_LP)');
        if isequal(cohortRoot, 0)
            error('extended_script11_dominant_period_run:NoRoot', 'No cohort folder selected.');
        end
    end

    paths = extended_script8_resolve_paths(cohortRoot);
    script11_assert_inputs_(paths);

    pal = extended_tol_bright_palette();
    theme = script11_theme_(cfg, pal);
    [outDirs, modeLabel] = script11_make_output_dirs_(cohortRoot, cfg.plotMode);

    logPath = fullfile(outDirs.logs, sprintf('Script11_DominantPeriod_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() script11_fclose_(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 11: Dominant UR periods (ClusterMembership) ===\n');
    fprintf('Cohort:  %s\n', paths.cohortRoot);
    fprintf('Output:  %s\n', outDirs.root);
    fprintf('Mode:    %s | dpi=%g | ext=%s\n', cfg.plotMode, theme.dpi, theme.ext);
    script11_log_(LOG, '%s v%s started', SCRIPT_NAME, SCRIPT_VERSION);
    script11_log_(LOG, 'Cohort: %s', paths.cohortRoot);
    script11_log_(LOG, 'profilesXlsx: %s', paths.profilesXlsx);
    script11_log_(LOG, 'Period source: Script 7 ClusterMembership.RawPeriod_h (no WP_TS)');

    clusterSummary = script11_read_sheet_(paths.profilesXlsx, 'ClusterSummary');
    clusterMembership = script11_read_sheet_(paths.profilesXlsx, 'ClusterMembership');
    if isempty(clusterSummary) || isempty(clusterMembership)
        error('extended_script11_dominant_period_run:MissingClusters', ...
            'ClusterSummary or ClusterMembership missing in %s', paths.profilesXlsx);
    end

    candTable = script11_prepare_membership_(clusterMembership, LOG);
    panelSpec = cfg.script11.panels;
    resolved = script11_resolve_panels_(clusterSummary, panelSpec, LOG);

    p0Ref = script11_resolve_p0_(paths.handoffRoot, cfg.script11.defaultP0_h, LOG);

    [clusterTable, mouseTable] = script11_compute_tables_(candTable, resolved, clusterSummary, p0Ref, cfg, LOG);

    settingsTable = script11_make_settings_table_(SCRIPT_NAME, SCRIPT_VERSION, paths, p0Ref, cfg);
    tablePaths = script11_write_tables_(outDirs.tables, settingsTable, clusterTable, mouseTable, LOG);

    figPath = fullfile(outDirs.figures, ['Supp_DominantPeriod_ClusterPeriod_3Cluster' theme.ext]);
    script11_plot_violin_panels_(figPath, mouseTable, candTable, resolved, theme, pal, LOG);
    script11_log_(LOG, 'Wrote figure %s', figPath);

    out = struct();
    out.cohortRoot = paths.cohortRoot;
    out.outputRoot = outDirs.root;
    out.modeLabel = modeLabel;
    out.clusterSummaryTable = clusterTable;
    out.perMouseTable = mouseTable;
    out.tablePaths = tablePaths;
    out.figurePath = figPath;
    out.logPath = logPath;
    out.p0Reference_h = p0Ref;

    fprintf('Cluster summary: %s\n', tablePaths.clusterSummaryCsv);
    fprintf('Per-mouse table: %s\n', tablePaths.perMouseCsv);
    fprintf('Workbook:        %s\n', tablePaths.xlsx);
    fprintf('Figure:          %s\n', figPath);
    fprintf('Log:             %s\n', logPath);
end

%% Table export
function settingsTable = script11_make_settings_table_(scriptName, scriptVersion, paths, p0Ref, cfg)
    names = ["ScriptName"; "ScriptVersion"; "RunTimestamp"; "CohortRoot"; ...
        "ProfilesXlsx"; "PlotMode"; "P0_Reference_h"; "PeriodSource"; "Claim"];
    vals = [string(scriptName); string(scriptVersion); string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')); ...
        string(paths.cohortRoot); string(paths.profilesXlsx); ...
        string(cfg.plotMode); string(p0Ref); ...
        "Script7 ClusterMembership.RawPeriod_h (candidate-level; no WP_TS)"; ...
        "Per-mouse mean+median of candidate periods within locked clusters; cohort mean+median across mice"];
    settingsTable = table(names, vals, 'VariableNames', {'Setting', 'Value'});
end

function tablePaths = script11_write_tables_(tablesDir, settingsTable, clusterTable, mouseTable, LOG)
    tablePaths = struct();
    tablePaths.clusterSummaryCsv = fullfile(tablesDir, 'DominantPeriod_ClusterSummary.csv');
    tablePaths.perMouseCsv = fullfile(tablesDir, 'DominantPeriod_PerMouse.csv');
    tablePaths.xlsx = fullfile(tablesDir, 'DominantPeriod_Output.xlsx');

    writetable(clusterTable, tablePaths.clusterSummaryCsv);
    writetable(mouseTable, tablePaths.perMouseCsv);
    script11_log_(LOG, 'Wrote %s (%d rows)', tablePaths.clusterSummaryCsv, height(clusterTable));
    script11_log_(LOG, 'Wrote %s (%d rows)', tablePaths.perMouseCsv, height(mouseTable));

    script11_safe_delete_file_(tablePaths.xlsx);
    script11_safe_writetable_(settingsTable, tablePaths.xlsx, 'Settings');
    script11_safe_writetable_(clusterTable, tablePaths.xlsx, 'ClusterSummary');
    script11_safe_writetable_(mouseTable, tablePaths.xlsx, 'PerMouse');
    script11_log_(LOG, 'Wrote %s (Settings, ClusterSummary, PerMouse)', tablePaths.xlsx);
end

function script11_safe_writetable_(T, xlsxPath, sheetName)
    if isempty(T)
        T = table();
    end
    sheetName = char(string(sheetName));
    if strlength(string(sheetName)) > 31
        sheetName = char(extractBefore(string(sheetName), 32));
    end
    try
        writetable(T, xlsxPath, 'Sheet', sheetName);
    catch ME
        warning('script11:WriteSheetFailed', 'Could not write sheet %s: %s', sheetName, ME.message);
    end
end

function script11_safe_delete_file_(filePath)
    if isfile(filePath)
        try
            delete(filePath);
        catch
        end
    end
end

%% Configuration / paths
function cfg = script11_fill_cfg_(cfg)
    if ~isfield(cfg, 'script11') || isempty(cfg.script11)
        cfg.script11 = struct();
    end
    if ~isfield(cfg.script11, 'panels') || isempty(cfg.script11.panels)
        cfg.script11.panels = [ ...
            struct('BandName', "UR_1_3", 'ClusterRank', 1, 'MainText', true); ...
            struct('BandName', "UR_3_6", 'ClusterRank', 1, 'MainText', true); ...
            struct('BandName', "UR_3_6", 'ClusterRank', 2, 'MainText', false)];
    end
    if ~isfield(cfg.script11, 'defaultP0_h') || isempty(cfg.script11.defaultP0_h)
        cfg.script11.defaultP0_h = 24.0;
    end
    if ~isfield(cfg.script11, 'harmonicKMax') || isempty(cfg.script11.harmonicKMax)
        cfg.script11.harmonicKMax = 12;
    end
end

function script11_assert_inputs_(paths)
    missing = strings(0, 1);
    if ~isfile(paths.profilesXlsx)
        missing(end + 1, 1) = "Script7 profiles xlsx"; %#ok<AGROW>
    end
    if ~isempty(missing)
        error('extended_script11_dominant_period_run:MissingInputs', ...
            'Missing: %s. Run Script 7 first.', strjoin(missing, ', '));
    end
end

function [outDirs, modeLabel] = script11_make_output_dirs_(cohortRoot, plotMode)
    switch lower(string(plotMode))
        case "publication"
            modeLabel = "Publication";
        otherwise
            modeLabel = "Development";
    end
    outDirs = struct();
    outDirs.root = fullfile(cohortRoot, sprintf('Script11_DominantPeriod_%s', char(modeLabel)));
    outDirs.tables = fullfile(outDirs.root, 'Tables');
    outDirs.figures = fullfile(outDirs.root, 'Figures');
    outDirs.logs = fullfile(outDirs.root, 'Logs');
    script11_ensure_dir_(outDirs.root);
    script11_ensure_dir_(outDirs.tables);
    script11_ensure_dir_(outDirs.figures);
    script11_ensure_dir_(outDirs.logs);
end

function theme = script11_theme_(cfg, pal)
    theme = struct();
    theme.dpi = cfg.plot.saveDpi;
    theme.ext = cfg.plot.figExt;
    theme.fontName = pal.fontName;
    theme.fontSize = 14;
    theme.labelSize = 16;
    theme.titleSize = 17;
    theme.axesLineWidth = pal.axesLineWidth;
    theme.tickDir = pal.tickDir;
    theme.medianColor = [0.85 0.00 0.55];  % magenta dotted — reference style
    theme.meanColor = [0.10 0.10 0.10];    % black solid — reference style
    theme.medianLineWidth = 1.5;
    theme.meanLineWidth = 1.6;
    theme.violinAlpha = 0.55;
    theme.violinWidth = 0.38;
end

%% Data load / membership
function T = script11_read_sheet_(xlsxPath, sheetName)
    T = table();
    if ~isfile(xlsxPath)
        return;
    end
    try
        opts = detectImportOptions(xlsxPath, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
        T = readtable(xlsxPath, opts);
    catch
        T = table();
    end
end

function cand = script11_prepare_membership_(clusterMembership, LOG)
    needed = {'CandidateID','SignalID','ClusterID','BandName','RawPeriod_h'};
    if ~all(ismember(needed, clusterMembership.Properties.VariableNames))
        error('extended_script11_dominant_period_run:BadMembership', ...
            'ClusterMembership missing required columns: %s', strjoin(needed, ', '));
    end
    cand = clusterMembership;
    cand.CandidateID = string(cand.CandidateID);
    cand.SignalID = string(cand.SignalID);
    cand.ClusterID = string(cand.ClusterID);
    cand.BandName = string(cand.BandName);
    cand.RawPeriod_h = double(cand.RawPeriod_h);
    if ismember('ClusterRank', cand.Properties.VariableNames)
        cand.ClusterRank = double(cand.ClusterRank);
    end
    if ismember('Photoperiod_h', cand.Properties.VariableNames)
        cand.Photoperiod_h = double(cand.Photoperiod_h);
    end
    keep = isfinite(cand.RawPeriod_h) & cand.RawPeriod_h > 0 & ...
        strlength(cand.SignalID) > 0 & strlength(cand.ClusterID) > 0;
    cand = cand(keep, :);
    if height(cand) == 0
        error('extended_script11_dominant_period_run:EmptyMembership', ...
            'No usable ClusterMembership rows with finite RawPeriod_h.');
    end
    script11_log_(LOG, 'ClusterMembership candidates: %d (mice=%d, clusters=%d)', ...
        height(cand), numel(unique(cand.SignalID)), numel(unique(cand.ClusterID)));
end

function resolved = script11_resolve_panels_(clusterSummary, panelSpec, LOG)
    resolved = struct('BandName', {}, 'ClusterRank', {}, 'MainText', {}, ...
        'ClusterID', {}, 'PeriodLow_h', {}, 'PeriodHigh_h', {}, ...
        'FilterLow_h', {}, 'FilterHigh_h', {}, 'PeriodCentre_h', {}, 'PanelIndex', {});

    CS = clusterSummary;
    CS.BandName = string(CS.BandName);
    CS.ClusterID = string(CS.ClusterID);

    for i = 1:numel(panelSpec)
        band = string(panelSpec(i).BandName);
        rank = double(panelSpec(i).ClusterRank);
        cid = script11_cluster_by_rank_(CS, band, rank);
        if strlength(cid) == 0
            error('extended_script11_dominant_period_run:ClusterMissing', ...
                'No cluster for band %s rank %d in ClusterSummary.', char(band), rank);
        end
        hit = CS(string(CS.ClusterID) == cid, :);
        if isempty(hit)
            error('extended_script11_dominant_period_run:ClusterRowMissing', ...
                'ClusterID %s not found in ClusterSummary.', char(cid));
        end
        script11_log_(LOG, 'Panel %d: %s rank %d -> %s', i, char(band), rank, char(cid));

        resolved(i).BandName = band;
        resolved(i).ClusterRank = rank;
        resolved(i).MainText = logical(panelSpec(i).MainText);
        resolved(i).ClusterID = cid;
        resolved(i).PeriodLow_h = double(hit.PeriodLow_h(1));
        resolved(i).PeriodHigh_h = double(hit.PeriodHigh_h(1));
        resolved(i).FilterLow_h = double(hit.FilterLow_h(1));
        resolved(i).FilterHigh_h = double(hit.FilterHigh_h(1));
        resolved(i).PeriodCentre_h = double(hit.PeriodCentre_h(1));
        resolved(i).PanelIndex = i;
    end
end

function cid = script11_cluster_by_rank_(CS, bandName, rank)
    cid = "";
    sub = CS(string(CS.BandName) == string(bandName), :);
    if isempty(sub), return; end
    rank = double(rank);
    if ismember('ClusterRank', sub.Properties.VariableNames)
        hit = sub(double(sub.ClusterRank) == rank, :);
        if ~isempty(hit)
            cid = string(hit.ClusterID(1));
            return;
        end
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

function p0 = script11_resolve_p0_(handoffRoot, defaultP0, LOG)
    p0 = defaultP0;
    if ~isfolder(handoffRoot)
        script11_log_(LOG, 'Handoff folder missing; using default P0=%.2f h', p0);
        return;
    end
    hits = dir(fullfile(handoffRoot, 'CoreSummary__*.mat'));
    if isempty(hits)
        script11_log_(LOG, 'No CoreSummary__*.mat; using default P0=%.2f h', p0);
        return;
    end
    try
        S = load(fullfile(hits(1).folder, hits(1).name));
        anchorVals = script11_collect_anchor_periods_(S);
        if ~isempty(anchorVals)
            p0 = median(anchorVals, 'omitnan');
            script11_log_(LOG, 'P0 reference from Core handoff median = %.3f h (n=%d files/mice aggregate)', ...
                p0, numel(anchorVals));
        else
            script11_log_(LOG, 'No anchor periods in handoff; using default P0=%.2f h', p0);
        end
    catch ME
        script11_log_(LOG, 'Could not load Core handoff for P0 (%s); default %.2f h', ME.message, p0);
    end
end

function vals = script11_collect_anchor_periods_(S)
    vals = [];
    if isfield(S, 'handoff') && isstruct(S.handoff)
        H = S.handoff;
        if isfield(H, 'hsub') && isstruct(H.hsub) && isfield(H.hsub, 'anchorPeriod_h')
            v = double(H.hsub.anchorPeriod_h(:));
            vals = [vals; v(isfinite(v) & v > 0)]; %#ok<AGROW>
        end
    end
    if isfield(S, 'anchorReport') && istable(S.anchorReport)
        AR = S.anchorReport;
        for col = {'P0_h', 'AnchorPeriod_h', 'P0w_median_h'}
            c = col{1};
            if ismember(c, AR.Properties.VariableNames)
                v = double(AR.(c));
                vals = [vals; v(isfinite(v) & v > 0)]; %#ok<AGROW>
            end
        end
    end
end

%% Statistics
function [clusterTable, mouseTable] = script11_compute_tables_(candTable, resolved, clusterSummary, p0Ref, cfg, LOG)
    mouseRows = {};
    clusterRows = {};

    CS = clusterSummary;
    CS.ClusterID = string(CS.ClusterID);

    for i = 1:numel(resolved)
        R = resolved(i);
        sub = candTable(string(candTable.ClusterID) == string(R.ClusterID), :);
        if height(sub) == 0
            error('extended_script11_dominant_period_run:NoCandidates', ...
                'No ClusterMembership candidates for cluster %s', char(R.ClusterID));
        end

        mice = unique(string(sub.SignalID));
        mice = mice(strlength(mice) > 0);
        perMouseMedian = NaN(numel(mice), 1);
        perMouseMean = NaN(numel(mice), 1);
        perMouseIQR = NaN(numel(mice), 1);
        perMouseSD = NaN(numel(mice), 1);
        perMouseCand = zeros(numel(mice), 1);

        for m = 1:numel(mice)
            ms = sub(string(sub.SignalID) == mice(m), :);
            y = double(ms.RawPeriod_h);
            y = y(isfinite(y) & y > 0);
            perMouseMedian(m) = median(y, 'omitnan');
            perMouseMean(m) = mean(y, 'omitnan');
            perMouseIQR(m) = script11_iqr_(y);
            if numel(y) > 1
                perMouseSD(m) = std(y);
            else
                perMouseSD(m) = 0;
            end
            perMouseCand(m) = numel(y);
            mouseRows(end + 1, :) = {string(R.ClusterID), string(R.BandName), double(R.ClusterRank), ...
                logical(R.MainText), string(mice(m)), perMouseMedian(m), perMouseMean(m), ...
                perMouseIQR(m), perMouseSD(m), double(perMouseCand(m))}; %#ok<AGROW>
        end

        medCohort = median(perMouseMedian, 'omitnan');
        meanCohort = mean(perMouseMean, 'omitnan');
        iqrMice = script11_iqr_(perMouseMedian);
        sdMice = std(perMouseMedian, 'omitnan');
        medWithin = median(perMouseIQR, 'omitnan');
        [nearHarm, deltaHarm, kNear] = script11_nearest_harmonic_(medCohort, p0Ref, cfg.script11.harmonicKMax);

        hit = CS(string(CS.ClusterID) == string(R.ClusterID), :);
        script7Median = NaN;
        script7IQR = NaN;
        if ~isempty(hit)
            if ismember('MedianRawPeriod_h', hit.Properties.VariableNames)
                script7Median = double(hit.MedianRawPeriod_h(1));
            end
            if ismember('IQRRawPeriod_h', hit.Properties.VariableNames)
                script7IQR = double(hit.IQRRawPeriod_h(1));
            end
        end

        clusterRows(end + 1, :) = {string(R.ClusterID), string(R.BandName), double(R.ClusterRank), ...
            logical(R.MainText), double(R.PeriodLow_h), double(R.PeriodHigh_h), ...
            double(R.PeriodCentre_h), medCohort, meanCohort, iqrMice, sdMice, medWithin, ...
            double(numel(mice)), double(height(unique(string(sub.CandidateID)))), ...
            double(sum(perMouseCand)), script7Median, script7IQR, ...
            p0Ref, kNear, nearHarm, deltaHarm}; %#ok<AGROW>

        script11_log_(LOG, '%s: median=%.3f h, mean=%.3f h, IQR mice=%.3f, n=%d mice', ...
            char(R.ClusterID), medCohort, meanCohort, iqrMice, numel(mice));
    end

    mouseHdr = {'ClusterID','BandName','ClusterRank','IncludeInMainText','SignalID', ...
        'MedianPeriod_h','MeanPeriod_h','IQR_Period_h','SD_Period_h','N_Candidates'};
    clusterHdr = {'ClusterID','BandName','ClusterRank','IncludeInMainText', ...
        'PeriodLow_h','PeriodHigh_h','PeriodCentre_h', ...
        'MedianTau_Cohort_h','MeanTau_Cohort_h','IQR_AcrossMice_h','SD_AcrossMice_h', ...
        'MedianWithinMouseIQR_h','N_Mice','N_Candidates','N_CandidatePeriodValues', ...
        'Script7_MedianRawPeriod_h','Script7_IQRRawPeriod_h', ...
        'P0_Reference_h','NearestHarmonic_k','NearestHarmonic_h','DeltaFromHarmonic_h'};

    mouseTable = script11_rows_to_table_(mouseRows, mouseHdr);
    clusterTable = script11_rows_to_table_(clusterRows, clusterHdr);
end

function [nearH, deltaH, kNear] = script11_nearest_harmonic_(period_h, p0, kMax)
    nearH = NaN;
    deltaH = NaN;
    kNear = NaN;
    if ~isfinite(period_h) || ~isfinite(p0) || p0 <= 0
        return;
    end
    ks = 2:double(kMax);
    harmonics = p0 ./ ks;
    [deltaH, idx] = min(abs(period_h - harmonics));
    nearH = harmonics(idx);
    kNear = ks(idx);
end

%% Plotting — per-mouse candidate-period violins (mean + median lines)
function script11_plot_violin_panels_(outPath, mouseTable, candTable, resolved, theme, pal, LOG)
    nPanels = numel(resolved);
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 300 * nPanels + 120 560]);
    t = tiledlayout(fig, 1, nPanels, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nPanels
        R = resolved(i);
        ax = nexttile(t, i);
        hold(ax, 'on'); set(ax, 'Color', 'w');

        subCand = candTable(string(candTable.ClusterID) == string(R.ClusterID), :);
        subMouse = mouseTable(string(mouseTable.ClusterID) == string(R.ClusterID), :);
        if isempty(subMouse)
            title(ax, sprintf('%s C%02d — no mice', char(R.BandName), R.ClusterRank), ...
                'Interpreter', 'none');
            continue;
        end

        [~, sortIdx] = sort(double(subMouse.MedianPeriod_h), 'ascend', 'MissingPlacement', 'last');
        subMouse = subMouse(sortIdx, :);
        mice = string(subMouse.SignalID);
        bandColor = script11_band_color_(pal, R.BandName);

        for m = 1:numel(mice)
            y = double(subCand.RawPeriod_h(string(subCand.SignalID) == mice(m)));
            y = y(isfinite(y) & y > 0);
            if isempty(y), continue; end
            script11_draw_violin_(ax, m, y, bandColor, theme);
        end

        % Cohort median / mean of per-mouse summaries (horizontal guides)
        medCohort = median(double(subMouse.MedianPeriod_h), 'omitnan');
        meanCohort = mean(double(subMouse.MeanPeriod_h), 'omitnan');
        if isfinite(medCohort)
            yline(ax, medCohort, ':', 'Color', theme.medianColor, ...
                'LineWidth', theme.medianLineWidth, 'HandleVisibility', 'off');
        end
        if isfinite(meanCohort)
            yline(ax, meanCohort, '-', 'Color', theme.meanColor, ...
                'LineWidth', theme.meanLineWidth, 'HandleVisibility', 'off');
        end

        ax.XLim = [0.4, max(numel(mice), 1) + 0.6];
        ax.XTick = 1:numel(mice);
        ax.XTickLabel = arrayfun(@(k) sprintf('M%d', k), 1:numel(mice), 'UniformOutput', false);
        xlabel(ax, 'Mouse (sorted by median τ)', 'FontName', theme.fontName, ...
            'FontSize', theme.labelSize, 'FontWeight', 'bold');

        yLo = R.FilterLow_h;
        yHi = R.FilterHigh_h;
        if ~isfinite(yLo) || ~isfinite(yHi) || yHi <= yLo
            yLo = R.PeriodLow_h;
            yHi = R.PeriodHigh_h;
        end
        pad = 0.08 * max(yHi - yLo, 0.1);
        ylim(ax, [yLo - pad, yHi + pad]);
        ylabel(ax, 'Candidate period (h)', 'FontName', theme.fontName, ...
            'FontSize', theme.labelSize, 'FontWeight', 'bold');
        title(ax, sprintf('%s C%02d (%.2f–%.2f h)', char(R.BandName), R.ClusterRank, ...
            R.PeriodLow_h, R.PeriodHigh_h), ...
            'FontName', theme.fontName, 'FontSize', theme.titleSize, 'FontWeight', 'bold', ...
            'Interpreter', 'none');
        script11_style_axes_(ax, theme);
    end

    exportgraphics(fig, outPath, 'Resolution', theme.dpi);
    close(fig);
    script11_log_(LOG, 'Figure exported: %s', outPath);
end

function script11_draw_violin_(ax, xCenter, y, faceColor, theme)
    y = y(:);
    y = y(isfinite(y));
    if numel(y) < 2
        if numel(y) == 1
            plot(ax, xCenter, y, 'o', 'Color', faceColor, 'MarkerFaceColor', faceColor, ...
                'MarkerSize', 5, 'HandleVisibility', 'off');
        end
        return;
    end

    nBins = max(8, min(24, round(sqrt(numel(y)))));
    yMin = min(y);
    yMax = max(y);
    if yMax <= yMin
        yMax = yMin + 0.01;
    end
    edges = linspace(yMin, yMax, nBins + 1);
    counts = histcounts(y, edges, 'Normalization', 'pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    maxC = max(counts);
    if ~isfinite(maxC) || maxC <= 0
        maxC = 1;
    end
    halfW = theme.violinWidth * (counts / maxC);
    xx = [xCenter - halfW, fliplr(xCenter + halfW * 0.05)];
    yy = [centers, fliplr(centers)];
    patch(ax, xx, yy, faceColor, 'FaceAlpha', theme.violinAlpha, ...
        'EdgeColor', faceColor * 0.65, 'LineWidth', 0.6, 'HandleVisibility', 'off');

    medY = median(y, 'omitnan');
    meanY = mean(y, 'omitnan');
    w = theme.violinWidth;
    xLine = [xCenter - w, xCenter + w];
    plot(ax, xLine, [medY, medY], ':', 'Color', theme.medianColor, ...
        'LineWidth', theme.medianLineWidth, 'HandleVisibility', 'off');
    plot(ax, xLine, [meanY, meanY], '-', 'Color', theme.meanColor, ...
        'LineWidth', theme.meanLineWidth, 'HandleVisibility', 'off');
end

function c = script11_band_color_(pal, bandName)
    bandName = char(string(bandName));
    if isKey(pal.band, bandName)
        c = pal.band(bandName);
    else
        c = pal.base(1, :);
    end
end

%% Utilities
function T = script11_rows_to_table_(rows, hdr)
    if isempty(rows)
        T = cell2table(cell(0, numel(hdr)), 'VariableNames', hdr);
        return;
    end
    T = cell2table(rows, 'VariableNames', hdr);
end

function q = script11_iqr_(x)
    x = x(isfinite(x));
    if isempty(x)
        q = NaN;
        return;
    end
    q = iqr(x);
end

function script11_ensure_dir_(d)
    if ~isfolder(d)
        mkdir(d);
    end
end

function script11_style_axes_(ax, theme)
    box(ax, 'off');
    set(ax, 'TickDir', theme.tickDir, 'FontName', theme.fontName, ...
        'FontSize', theme.fontSize, 'LineWidth', theme.axesLineWidth);
end

function script11_log_(LOG, fmt, varargin)
    line = sprintf(fmt, varargin{:});
    fprintf('%s\n', line);
    if ~isempty(LOG) && LOG > 0
        fprintf(LOG, '%s\n', line);
    end
end

function script11_fclose_(LOG)
    if ~isempty(LOG) && LOG > 0
        fclose(LOG);
    end
end
