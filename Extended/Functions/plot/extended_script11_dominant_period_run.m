function out = extended_script11_dominant_period_run(cohortRoot, cfg)
%EXTENDED_SCRIPT11_DOMINANT_PERIOD_RUN Dominant cluster period tables + supp boxplots.
%
%   Inputs (read-only; no WP_TS):
%     Script 7  SelectedValidatedUR_*Profiles/Tables/*_Output.xlsx
%               ClusterSummary + ClusterMembership (RawPeriod_h per candidate)
%
%   Per mouse: mean + median of candidate RawPeriod_h within each locked cluster.
%   Cohort: mean + median across mice (of those per-mouse summaries).
%   Sex inferred from SignalID (same rules as Script 12); figures coloured F/M.
%   Figures: KDE violins only when n≥6; n≤5 = points + mean/median ticks (no invented density).
%
%   Outputs under {cohortRoot}/Script11_DominantPeriod_{Publication|Development}/
%     Tables/DominantPeriod_ClusterSummary.csv
%     Tables/DominantPeriod_PerMouse.csv
%     Tables/DominantPeriod_ByCluster.csv   (Cluster, N, Median_h, Median_SD, Mean_h, Mean_SD)
%     Tables/DominantPeriod_ByMouse.csv    (Cluster, BandName, Median_h, Median_SD, Mean_h, Mean_SD, SignalID, Sex)
%     Tables/DominantPeriod_Output.xlsx
%     Figures/Supp_DominantPeriod_ClusterPeriod_{Band}_C{rank}.*   (stable mouse order; F/M coloured)
%     Figures/Supp_DominantPeriod_PopulationByCluster.*            (pooled; one box per cluster)
%     Figures/Supp_DominantPeriod_PopulationByCluster_BySex.*      (F|M boxes per cluster)

    SCRIPT_NAME = 'extended_script11_dominant_period_run';
    SCRIPT_VERSION = '2.4';

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
    activityZT = script11_read_sheet_(paths.profilesXlsx, 'ActivityComponent_24h');
    [clusterCompact, mouseCompact] = script11_make_compact_tables_(mouseTable, resolved, activityZT, LOG);

    settingsTable = script11_make_settings_table_(SCRIPT_NAME, SCRIPT_VERSION, paths, p0Ref, cfg);
    tablePaths = script11_write_tables_(outDirs.tables, settingsTable, clusterTable, mouseTable, ...
        clusterCompact, mouseCompact, LOG);

    figPaths = script11_write_all_figures_(outDirs.figures, mouseTable, candTable, resolved, theme, pal, LOG);

    out = struct();
    out.cohortRoot = paths.cohortRoot;
    out.outputRoot = outDirs.root;
    out.modeLabel = modeLabel;
    out.clusterSummaryTable = clusterTable;
    out.perMouseTable = mouseTable;
    out.byClusterTable = clusterCompact;
    out.byMouseTable = mouseCompact;
    out.tablePaths = tablePaths;
    out.figurePaths = figPaths;
    out.logPath = logPath;
    out.p0Reference_h = p0Ref;

    fprintf('Cluster summary: %s\n', tablePaths.clusterSummaryCsv);
    fprintf('Per-mouse table: %s\n', tablePaths.perMouseCsv);
    fprintf('By-cluster:      %s\n', tablePaths.byClusterCsv);
    fprintf('By-mouse:        %s\n', tablePaths.byMouseCsv);
    fprintf('Workbook:        %s\n', tablePaths.xlsx);
    for fi = 1:numel(figPaths)
        fprintf('Figure:          %s\n', figPaths(fi));
    end
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
        "ByCluster.N=mice with ridge candidates; N_Mice_Activity from ActivityComponent (may be lower). ByMouse.N_Candidates=period values (duplicates overlap on plot). Activity panel n ≠ Script11 period N."];
    settingsTable = table(names, vals, 'VariableNames', {'Setting', 'Value'});
end

function tablePaths = script11_write_tables_(tablesDir, settingsTable, clusterTable, mouseTable, ...
        clusterCompact, mouseCompact, LOG)
    tablePaths = struct();
    tablePaths.clusterSummaryCsv = fullfile(tablesDir, 'DominantPeriod_ClusterSummary.csv');
    tablePaths.perMouseCsv = fullfile(tablesDir, 'DominantPeriod_PerMouse.csv');
    tablePaths.byClusterCsv = fullfile(tablesDir, 'DominantPeriod_ByCluster.csv');
    tablePaths.byMouseCsv = fullfile(tablesDir, 'DominantPeriod_ByMouse.csv');
    tablePaths.xlsx = fullfile(tablesDir, 'DominantPeriod_Output.xlsx');

    writetable(clusterTable, tablePaths.clusterSummaryCsv);
    writetable(mouseTable, tablePaths.perMouseCsv);
    writetable(clusterCompact, tablePaths.byClusterCsv);
    writetable(mouseCompact, tablePaths.byMouseCsv);
    script11_log_(LOG, 'Wrote %s (%d rows)', tablePaths.clusterSummaryCsv, height(clusterTable));
    script11_log_(LOG, 'Wrote %s (%d rows)', tablePaths.perMouseCsv, height(mouseTable));
    script11_log_(LOG, 'Wrote %s (%d rows)', tablePaths.byClusterCsv, height(clusterCompact));
    script11_log_(LOG, 'Wrote %s (%d rows)', tablePaths.byMouseCsv, height(mouseCompact));

    script11_safe_delete_file_(tablePaths.xlsx);
    script11_safe_writetable_(settingsTable, tablePaths.xlsx, 'Settings');
    script11_safe_writetable_(clusterTable, tablePaths.xlsx, 'ClusterSummary');
    script11_safe_writetable_(mouseTable, tablePaths.xlsx, 'PerMouse');
    script11_safe_writetable_(clusterCompact, tablePaths.xlsx, 'ByCluster');
    script11_safe_writetable_(mouseCompact, tablePaths.xlsx, 'ByMouse');
    script11_log_(LOG, 'Wrote %s (Settings, ClusterSummary, PerMouse, ByCluster, ByMouse)', tablePaths.xlsx);
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
    theme.violinNGrid = 401;               % KDE evaluation points
    theme.violinBwPad = 3.0;               % support ≈ this × bw so density → 0 (pointed tips)
    theme.violinBwCapFrac = 0.50;          % mild bw cap vs span (smooth, not wild)
    theme.violinMinN = 6;                  % violins only for n≥6; else points + ticks
    theme.showCandidatePoints = true;      % jitter raw periods on mouse violins
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
        sexVec = script11_infer_sex_(mice);
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
                logical(R.MainText), string(mice(m)), string(sexVec(m)), perMouseMedian(m), perMouseMean(m), ...
                perMouseIQR(m), perMouseSD(m), double(perMouseCand(m))}; %#ok<AGROW>
        end

        medCohort = median(perMouseMedian, 'omitnan');
        meanCohort = mean(perMouseMean, 'omitnan');
        iqrMice = script11_iqr_(perMouseMedian);
        sdMice = std(perMouseMedian, 'omitnan');
        medWithin = median(perMouseIQR, 'omitnan');
        [nearHarm, deltaHarm, kNear] = script11_nearest_harmonic_(medCohort, p0Ref, cfg.script11.harmonicKMax);

        nF = sum(sexVec == "Female");
        nM = sum(sexVec == "Male");
        medF = median(perMouseMedian(sexVec == "Female"), 'omitnan');
        medM = median(perMouseMedian(sexVec == "Male"), 'omitnan');
        meanF = mean(perMouseMean(sexVec == "Female"), 'omitnan');
        meanM = mean(perMouseMean(sexVec == "Male"), 'omitnan');

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
            double(numel(mice)), double(nF), double(nM), medF, medM, meanF, meanM, ...
            double(height(unique(string(sub.CandidateID)))), ...
            double(sum(perMouseCand)), script7Median, script7IQR, ...
            p0Ref, kNear, nearHarm, deltaHarm}; %#ok<AGROW>

        script11_log_(LOG, '%s: median=%.3f h, mean=%.3f h, IQR mice=%.3f, n=%d (F=%d M=%d)', ...
            char(R.ClusterID), medCohort, meanCohort, iqrMice, numel(mice), nF, nM);
    end

    mouseHdr = {'ClusterID','BandName','ClusterRank','IncludeInMainText','SignalID','Sex', ...
        'MedianPeriod_h','MeanPeriod_h','IQR_Period_h','SD_Period_h','N_Candidates'};
    clusterHdr = {'ClusterID','BandName','ClusterRank','IncludeInMainText', ...
        'PeriodLow_h','PeriodHigh_h','PeriodCentre_h', ...
        'MedianTau_Cohort_h','MeanTau_Cohort_h','IQR_AcrossMice_h','SD_AcrossMice_h', ...
        'MedianWithinMouseIQR_h','N_Mice','N_Female','N_Male', ...
        'MedianTau_Female_h','MedianTau_Male_h','MeanTau_Female_h','MeanTau_Male_h', ...
        'N_Candidates','N_CandidatePeriodValues', ...
        'Script7_MedianRawPeriod_h','Script7_IQRRawPeriod_h', ...
        'P0_Reference_h','NearestHarmonic_k','NearestHarmonic_h','DeltaFromHarmonic_h'};

    mouseTable = script11_rows_to_table_(mouseRows, mouseHdr);
    clusterTable = script11_rows_to_table_(clusterRows, clusterHdr);
end

function [clusterCompact, mouseCompact] = script11_make_compact_tables_(mouseTable, resolved, activityZT, LOG)
% Manuscript-oriented exports.
% ByCluster: Median_h/Median_SD = median and SD of per-mouse medians;
%            Mean_h/Mean_SD = mean and SD of per-mouse means; N = mice with period candidates.
%            N_Mice_Activity = unique mice in ActivityComponent_24h (may differ; see Settings).
% ByMouse:   Median_h/Mean_h from candidate periods; Median_SD and Mean_SD are both
%            the within-mouse SD of candidate periods; N_Candidates = plotted period count.

    if nargin < 3
        activityZT = table();
    end
    if nargin < 4
        LOG = [];
    end

    clusterRows = {};
    for i = 1:numel(resolved)
        R = resolved(i);
        sub = mouseTable(string(mouseTable.ClusterID) == string(R.ClusterID), :);
        meds = double(sub.MedianPeriod_h);
        means = double(sub.MeanPeriod_h);
        meds = meds(isfinite(meds));
        means = means(isfinite(means));
        n = height(sub);
        nCand = 0;
        if ismember('N_Candidates', sub.Properties.VariableNames)
            nCand = sum(double(sub.N_Candidates), 'omitnan');
        end
        nAct = script11_n_mice_activity_(activityZT, R.ClusterID);
        medH = NaN; medSD = NaN; meanH = NaN; meanSD = NaN;
        if ~isempty(meds)
            medH = median(meds, 'omitnan');
            if numel(meds) > 1
                medSD = std(meds, 'omitnan');
            else
                medSD = 0;
            end
        end
        if ~isempty(means)
            meanH = mean(means, 'omitnan');
            if numel(means) > 1
                meanSD = std(means, 'omitnan');
            else
                meanSD = 0;
            end
        end
        clusterRows(end + 1, :) = {string(script11_cluster_axis_label_(R)), double(n), ...
            medH, medSD, meanH, meanSD, double(nCand), double(nAct)}; %#ok<AGROW>
        if ~isempty(LOG)
            script11_log_(LOG, 'ByCluster %s: N_Mice=%d, N_Candidates=%d, N_Mice_Activity=%d', ...
                char(script11_cluster_axis_label_(R)), n, nCand, nAct);
        end
    end
    clusterCompact = script11_rows_to_table_(clusterRows, ...
        {'Cluster','N','Median_h','Median_SD','Mean_h','Mean_SD','N_Candidates','N_Mice_Activity'});

    mouseRows = {};
    if ~isempty(mouseTable)
        for r = 1:height(mouseTable)
            Rfake = struct('BandName', mouseTable.BandName(r), 'ClusterRank', mouseTable.ClusterRank(r));
            clusterLbl = string(script11_cluster_axis_label_(Rfake));
            sdWithin = double(mouseTable.SD_Period_h(r));
            nCand = NaN;
            if ismember('N_Candidates', mouseTable.Properties.VariableNames)
                nCand = double(mouseTable.N_Candidates(r));
            end
            mouseRows(end + 1, :) = {clusterLbl, string(mouseTable.BandName(r)), ...
                double(mouseTable.MedianPeriod_h(r)), sdWithin, ...
                double(mouseTable.MeanPeriod_h(r)), sdWithin, ...
                string(mouseTable.SignalID(r)), string(mouseTable.Sex(r)), nCand}; %#ok<AGROW>
        end
    end
    mouseCompact = script11_rows_to_table_(mouseRows, ...
        {'Cluster','BandName','Median_h','Median_SD','Mean_h','Mean_SD','SignalID','Sex','N_Candidates'});
end

function n = script11_n_mice_activity_(activityZT, clusterID)
    n = NaN;
    if isempty(activityZT)
        return;
    end
    need = {'ClusterID','SignalID'};
    if ~all(ismember(need, activityZT.Properties.VariableNames))
        return;
    end
    B = activityZT(string(activityZT.ClusterID) == string(clusterID), :);
    if isempty(B)
        n = 0;
        return;
    end
    ids = string(B.SignalID);
    ids = ids(strlength(ids) > 0);
    n = numel(unique(ids));
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

%% Plotting — separate per-cluster mouse violins + population-by-cluster
function figPaths = script11_write_all_figures_(figDir, mouseTable, candTable, resolved, theme, pal, LOG)
    figPaths = strings(0, 1);
    for i = 1:numel(resolved)
        R = resolved(i);
        stem = sprintf('Supp_DominantPeriod_ClusterPeriod_%s_C%02d', ...
            char(regexprep(char(R.BandName), '[^A-Za-z0-9]+', '_')), R.ClusterRank);
        outPath = fullfile(figDir, [stem theme.ext]);
        script11_plot_one_cluster_mice_(outPath, mouseTable, candTable, R, theme, pal);
        script11_log_(LOG, 'Wrote figure %s', outPath);
        figPaths(end + 1, 1) = string(outPath); %#ok<AGROW>
    end

    popPooled = fullfile(figDir, ['Supp_DominantPeriod_PopulationByCluster' theme.ext]);
    script11_plot_population_pooled_(popPooled, mouseTable, resolved, theme, pal);
    script11_log_(LOG, 'Wrote figure %s', popPooled);
    figPaths(end + 1, 1) = string(popPooled);

    popSex = fullfile(figDir, ['Supp_DominantPeriod_PopulationByCluster_BySex' theme.ext]);
    script11_plot_population_by_sex_(popSex, mouseTable, resolved, theme, pal);
    script11_log_(LOG, 'Wrote figure %s', popSex);
    figPaths(end + 1, 1) = string(popSex);
end

function script11_plot_one_cluster_mice_(outPath, mouseTable, candTable, R, theme, pal)
% Stable SignalID order (not sorted by sex or median); violins coloured by sex.
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 960 640]);
    ax = axes(fig); hold(ax, 'on'); set(ax, 'Color', 'w');

    subCand = candTable(string(candTable.ClusterID) == string(R.ClusterID), :);
    subMouse = mouseTable(string(mouseTable.ClusterID) == string(R.ClusterID), :);
    if isempty(subMouse)
        title(ax, sprintf('%s — no mice', script11_cluster_axis_label_(R)), 'Interpreter', 'none');
        exportgraphics(fig, outPath, 'Resolution', theme.dpi);
        close(fig);
        return;
    end

    if ~ismember('Sex', subMouse.Properties.VariableNames)
        subMouse.Sex = script11_infer_sex_(string(subMouse.SignalID));
    end

    % Stable alphabetical SignalID order so M# is reproducible across runs
    [~, sortIdx] = sort(string(subMouse.SignalID));
    subMouse = subMouse(sortIdx, :);
    mice = string(subMouse.SignalID);

    ySeen = [];
    for m = 1:numel(mice)
        y = double(subCand.RawPeriod_h(string(subCand.SignalID) == mice(m)));
        y = y(isfinite(y) & y > 0);
        if isempty(y), continue; end
        faceColor = script11_sex_color_(pal, subMouse.Sex(m));
        yExt = script11_draw_violin_(ax, m, y, faceColor, theme);
        ySeen = [ySeen; y(:); yExt(:)]; %#ok<AGROW>
    end

    ax.XLim = [0.4, max(numel(mice), 1) + 0.6];
    ax.XTick = 1:numel(mice);
    ax.XTickLabel = arrayfun(@(k) sprintf('M%d', k), 1:numel(mice), 'UniformOutput', false);
    xlabel(ax, 'Mouse', 'FontName', theme.fontName, ...
        'FontSize', theme.labelSize, 'FontWeight', 'bold');

    yHi = R.FilterHigh_h;
    if ~isfinite(yHi)
        yHi = R.PeriodHigh_h;
    end
    if ~isempty(ySeen)
        yHi = max(yHi, max(ySeen, [], 'omitnan'));
    end
    if ~isfinite(yHi) || yHi <= 0
        yHi = 1;
    end
    pad = 0.06 * max(yHi, 0.1);
    ylim(ax, [0, yHi + pad]);
    ylabel(ax, 'Period (h)', 'FontName', theme.fontName, ...
        'FontSize', theme.labelSize, 'FontWeight', 'bold');
    title(ax, sprintf('%s (%.2f–%.2f h)', script11_cluster_axis_label_(R), ...
        R.PeriodLow_h, R.PeriodHigh_h), ...
        'FontName', theme.fontName, 'FontSize', theme.titleSize, 'FontWeight', 'bold', ...
        'Interpreter', 'none');
    script11_style_axes_(ax, theme);
    script11_add_sex_mean_median_legend_(ax, theme, pal);

    exportgraphics(fig, outPath, 'Resolution', theme.dpi);
    close(fig);
end

function script11_plot_population_pooled_(outPath, mouseTable, resolved, theme, pal)
% One violin per cluster = all mice (no sex split).
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 780 640]);
    ax = axes(fig); hold(ax, 'on'); set(ax, 'Color', 'w');

    nC = numel(resolved);
    xLabels = strings(nC, 1);
    yAll = [];
    themePop = theme;
    themePop.showCandidatePoints = false;

    for i = 1:nC
        R = resolved(i);
        xLabels(i) = string(script11_cluster_axis_label_(R));
        subMouse = mouseTable(string(mouseTable.ClusterID) == string(R.ClusterID), :);
        y = double(subMouse.MedianPeriod_h);
        y = y(isfinite(y) & y > 0);
        if isempty(y)
            continue;
        end
        bandColor = extended_cluster_colour(pal, 'BandName', R.BandName, 'ClusterRank', R.ClusterRank);
        yExt = script11_draw_violin_(ax, i, y, bandColor, themePop);
        rng(11 + i);
        jitter = (rand(numel(y), 1) - 0.5) * 0.18;
        scatter(ax, i + jitter, y, 28, ...
            'MarkerFaceColor', bandColor, 'MarkerEdgeColor', bandColor * 0.7, ...
            'MarkerFaceAlpha', 0.55, 'HandleVisibility', 'off');
        yAll = [yAll; y(:); yExt(:)]; %#ok<AGROW>
    end

    set(ax, 'XTick', 1:nC, 'XTickLabel', cellstr(xLabels), 'XLim', [0.4, nC + 0.6]);
    xlabel(ax, 'Cluster', 'FontName', theme.fontName, 'FontSize', theme.labelSize, 'FontWeight', 'bold');
    ylabel(ax, 'Period (h)', 'FontName', theme.fontName, ...
        'FontSize', theme.labelSize, 'FontWeight', 'bold');
    title(ax, 'Population spread of per-mouse median period by cluster', ...
        'FontName', theme.fontName, 'FontSize', theme.titleSize, 'FontWeight', 'bold', ...
        'Interpreter', 'none');
    if ~isempty(yAll)
        yHi = max(yAll, [], 'omitnan');
        pad = 0.08 * max(yHi, 0.2);
        ylim(ax, [0, yHi + pad]);
    else
        ylim(ax, [0, 1]);
    end
    script11_style_axes_(ax, theme);
    script11_add_mean_median_legend_(ax, theme);

    exportgraphics(fig, outPath, 'Resolution', theme.dpi);
    close(fig);
end

function script11_plot_population_by_sex_(outPath, mouseTable, resolved, theme, pal)
% Per cluster: Female | Male violins of per-mouse median periods (exploratory).
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 920 640]);
    ax = axes(fig); hold(ax, 'on'); set(ax, 'Color', 'w');

    nC = numel(resolved);
    xLabels = strings(nC, 1);
    yAll = [];
    xLayout = extended_grouped_x_layout(2, 'HalfWidth', 0.22, 'Gap', 0.10);
    sexes = ["Female", "Male"];
    offsets = xLayout.offsets;

    if ~ismember('Sex', mouseTable.Properties.VariableNames)
        mouseTable.Sex = script11_infer_sex_(string(mouseTable.SignalID));
    end

    themePop = theme;
    themePop.showCandidatePoints = false;
    themePop.violinWidth = min(theme.violinWidth, xLayout.halfWidth * 2);

    for i = 1:nC
        R = resolved(i);
        xLabels(i) = string(script11_cluster_axis_label_(R));
        subMouse = mouseTable(string(mouseTable.ClusterID) == string(R.ClusterID), :);
        for si = 1:numel(sexes)
            y = double(subMouse.MedianPeriod_h(string(subMouse.Sex) == sexes(si)));
            y = y(isfinite(y) & y > 0);
            if isempty(y)
                continue;
            end
            x = i + offsets(si);
            faceColor = script11_sex_color_(pal, sexes(si));
            yExt = script11_draw_violin_(ax, x, y, faceColor, themePop);
            rng(11 + 10 * i + si);
            jitter = (rand(numel(y), 1) - 0.5) * 0.08;
            scatter(ax, x + jitter, y, 28, ...
                'MarkerFaceColor', faceColor, 'MarkerEdgeColor', faceColor * 0.65, ...
                'MarkerFaceAlpha', 0.65, 'HandleVisibility', 'off');
            yAll = [yAll; y(:); yExt(:)]; %#ok<AGROW>
        end
    end

    set(ax, 'XTick', 1:nC, 'XTickLabel', cellstr(xLabels), 'XLim', [0.35, nC + 0.65]);
    xlabel(ax, 'Cluster', 'FontName', theme.fontName, ...
        'FontSize', theme.labelSize, 'FontWeight', 'bold');
    ylabel(ax, 'Candidate period (h)', 'FontName', theme.fontName, ...
        'FontSize', theme.labelSize, 'FontWeight', 'bold');
    title(ax, 'Population spread of per-mouse median period by cluster and sex', ...
        'FontName', theme.fontName, 'FontSize', theme.titleSize, 'FontWeight', 'bold', ...
        'Interpreter', 'none');
    if ~isempty(yAll)
        yHi = max(yAll, [], 'omitnan');
        pad = 0.08 * max(yHi, 0.2);
        ylim(ax, [0, yHi + pad]);
    else
        ylim(ax, [0, 1]);
    end
    script11_style_axes_(ax, theme);
    script11_add_sex_mean_median_legend_(ax, theme, pal);

    exportgraphics(fig, outPath, 'Resolution', theme.dpi);
    close(fig);
end

function lbl = script11_cluster_axis_label_(R)
% e.g. UR_1_3 C01 -> "UR 1-3 C01"
    band = char(string(R.BandName));
    band = regexprep(band, '^UR_', 'UR ');
    band = strrep(band, '_', '-');
    lbl = sprintf('%s C%02d', band, double(R.ClusterRank));
end

function script11_add_mean_median_legend_(ax, theme)
    hMed = plot(ax, nan, nan, ':', 'Color', theme.medianColor, ...
        'LineWidth', theme.medianLineWidth, 'DisplayName', 'Median');
    hMean = plot(ax, nan, nan, '-', 'Color', theme.meanColor, ...
        'LineWidth', theme.meanLineWidth, 'DisplayName', 'Mean');
    lg = legend(ax, [hMed, hMean], {'Median', 'Mean'}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal');
    lg.Box = 'off';
    lg.FontName = theme.fontName;
    lg.FontSize = theme.fontSize;
end

function script11_add_sex_mean_median_legend_(ax, theme, pal)
    hF = patch(ax, nan, nan, pal.female, 'FaceAlpha', theme.violinAlpha, ...
        'EdgeColor', 'none', 'DisplayName', 'Female');
    hM = patch(ax, nan, nan, pal.male, 'FaceAlpha', theme.violinAlpha, ...
        'EdgeColor', 'none', 'DisplayName', 'Male');
    hMed = plot(ax, nan, nan, ':', 'Color', theme.medianColor, ...
        'LineWidth', theme.medianLineWidth, 'DisplayName', 'Median');
    hMean = plot(ax, nan, nan, '-', 'Color', theme.meanColor, ...
        'LineWidth', theme.meanLineWidth, 'DisplayName', 'Mean');
    lg = legend(ax, [hF, hM, hMed, hMean], {'Female', 'Male', 'Median', 'Mean'}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 4);
    lg.Box = 'off';
    lg.FontName = theme.fontName;
    lg.FontSize = theme.fontSize;
end

function yExt = script11_draw_violin_(ax, xCenter, y, faceColor, theme)
% Smooth KDE violin when n≥6 and periods are not identical; else points + mean/median ticks.
% Support extends ~3× bandwidth past the data so tips taper to a point (not flat trim).
% Bandwidth mildly capped vs span to limit empty whiskers; ylim may include short tips.
    yExt = [NaN; NaN];
    y = y(:);
    y = y(isfinite(y));
    if isempty(y)
        return;
    end

    medY = median(y, 'omitnan');
    meanY = mean(y, 'omitnan');
    w = theme.violinWidth * 0.85;
    xLine = [xCenter - w, xCenter + w];
    yMin = min(y);
    yMax = max(y);
    spanRaw = yMax - yMin;
    yExt = [yMin; yMax];

    minN = 6;
    if isfield(theme, 'violinMinN') && isfinite(theme.violinMinN)
        minN = max(3, round(double(theme.violinMinN)));
    end
    drawPointsOnly = numel(y) < minN || spanRaw < 1e-6;
    if drawPointsOnly
        rng(abs(round(1e3 * xCenter + sum(y))));
        if numel(y) == 1
            plot(ax, xCenter, y, 'o', 'Color', faceColor, 'MarkerFaceColor', faceColor, ...
                'MarkerSize', 6, 'HandleVisibility', 'off');
        else
            jitter = (rand(numel(y), 1) - 0.5) * theme.violinWidth * 0.55;
            scatter(ax, xCenter + jitter, y, 18, ...
                'MarkerFaceColor', faceColor, 'MarkerEdgeColor', faceColor * 0.55, ...
                'MarkerFaceAlpha', 0.75, 'HandleVisibility', 'off');
        end
        plot(ax, xLine, [medY, medY], ':', 'Color', theme.medianColor, ...
            'LineWidth', theme.medianLineWidth, 'HandleVisibility', 'off');
        plot(ax, xLine, [meanY, meanY], '-', 'Color', theme.meanColor, ...
            'LineWidth', theme.meanLineWidth, 'HandleVisibility', 'off');
        return;
    end

    bw = NaN;
    try
        [~, ~, bw] = ksdensity(y);
    catch
        bw = NaN;
    end
    if ~isfinite(bw) || bw <= 0
        bw = max(0.08 * spanRaw, 0.03);
    end
    bwCapFrac = 0.50;
    if isfield(theme, 'violinBwCapFrac') && isfinite(theme.violinBwCapFrac)
        bwCapFrac = double(theme.violinBwCapFrac);
    end
    bwCap = max(bwCapFrac * spanRaw, 0.05);
    bw = min(bw, bwCap);

    bwPad = 3.0;
    if isfield(theme, 'violinBwPad') && isfinite(theme.violinBwPad)
        bwPad = double(theme.violinBwPad);
    end
    % Enough margin that density ≈ 0 at the ends → pointed (non-flat) tips
    margin = max(bwPad * bw, 0.20 * spanRaw);
    nGrid = max(201, round(double(theme.violinNGrid)));
    yGrid = linspace(yMin - margin, yMax + margin, nGrid);

    f = [];
    yg = yGrid;
    try
        [f, yg] = ksdensity(y, yGrid, 'Bandwidth', bw);
    catch
        f = [];
    end
    if isempty(f) || ~any(isfinite(f) & f > 0)
        nBins = max(12, min(36, round(sqrt(numel(y)) * 3)));
        edges = linspace(yMin - margin, yMax + margin, nBins + 1);
        counts = histcounts(y, edges, 'Normalization', 'pdf');
        yg = (edges(1:end-1) + edges(2:end)) / 2;
        f = counts;
    end
    f = max(double(f(:))', 0);
    yg = double(yg(:))';
    maxF = max(f);
    if ~isfinite(maxF) || maxF <= 0
        maxF = 1;
    end
    f(1) = 0;
    f(end) = 0;
    halfW = theme.violinWidth * (f / maxF);
    xx = [xCenter - halfW, fliplr(xCenter + halfW)];
    yy = [yg, fliplr(yg)];
    patch(ax, xx, yy, faceColor, 'FaceAlpha', theme.violinAlpha, ...
        'EdgeColor', faceColor * 0.65, 'LineWidth', 0.6, 'HandleVisibility', 'off');
    % Include short tips in extent so ylim does not clip pointed ends
    yExt = [min(yg); max(yg)];

    if isfield(theme, 'showCandidatePoints') && theme.showCandidatePoints
        rng(abs(round(1e3 * xCenter + sum(y))));
        jitter = (rand(numel(y), 1) - 0.5) * theme.violinWidth * 0.55;
        scatter(ax, xCenter + jitter, y, 12, ...
            'MarkerFaceColor', faceColor * 0.55 + [1 1 1] * 0.45, ...
            'MarkerEdgeColor', faceColor * 0.5, 'MarkerFaceAlpha', 0.45, ...
            'HandleVisibility', 'off');
    end

    plot(ax, xLine, [medY, medY], ':', 'Color', theme.medianColor, ...
        'LineWidth', theme.medianLineWidth, 'HandleVisibility', 'off');
    plot(ax, xLine, [meanY, meanY], '-', 'Color', theme.meanColor, ...
        'LineWidth', theme.meanLineWidth, 'HandleVisibility', 'off');
end

function c = script11_sex_color_(pal, sex)
    sx = string(sex);
    if sx == "Female"
        c = pal.female;
    elseif sx == "Male"
        c = pal.male;
    else
        c = pal.pooled;
    end
end

function sx = script11_infer_sex_(signalID)
% Same SignalID heuristics as Script 12 (exploratory display / table annotation).
    s = lower(string(signalID));
    sx = repmat("Unknown", numel(s), 1);
    isFemale = contains(s, "female") | contains(s, "_f_") | contains(s, "-f-") | ...
        startsWith(s, "f_") | startsWith(s, "f-") | endsWith(s, "_f") | endsWith(s, "-f");
    isMale = contains(s, "male") | contains(s, "_m_") | contains(s, "-m-") | ...
        startsWith(s, "m_") | startsWith(s, "m-") | endsWith(s, "_m") | endsWith(s, "-m");
    isMale = isMale & ~contains(s, "female");
    sx(isMale) = "Male";
    sx(isFemale) = "Female";
end

function c = script11_band_color_(pal, bandName)
% Legacy wrapper — prefer extended_cluster_colour with ClusterRank when known.
    c = extended_cluster_colour(pal, 'BandName', bandName, 'ClusterRank', 1);
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
