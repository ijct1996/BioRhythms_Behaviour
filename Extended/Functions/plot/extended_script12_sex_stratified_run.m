function out = extended_script12_sex_stratified_run(cohortRoot, cfg)
%EXTENDED_SCRIPT12_SEX_STRATIFIED_RUN Sex-stratified activity grids + pre/post-LD amplitude.
%
%   For each locked cluster (UR_1_3 C01, UR_3_6 C01, UR_3_6 C02):
%     - Activity: Female | Male side-by-side L12–L22 (2×3 each)
%     - Amplitude: F/M overlaid vs photoperiod — pre-LD and post-LD first local max
%       minus baseline (search window ≈ 1× cluster period, capped)
%     - F vs M Mann–Whitney U per photoperiod; BH-FDR within ClusterID×Metric;
%       stars from Q_BH (* / ** / ***)
%
%   No phase coherence / ridge-power panels. Inputs: Script 7 profiles only.

    SCRIPT_NAME = 'extended_script12_sex_stratified_run';
    SCRIPT_VERSION = '1.7';

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    if ~isfield(cfg, 'plotMode') || isempty(cfg.plotMode)
        cfg.plotMode = 'development';
    end
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);
    cfg = script12_fill_cfg_(cfg);

    if nargin < 1 || isempty(cohortRoot)
        cohortRoot = uigetdir(pwd, 'Select cohort results folder (e.g. C57_LP)');
        if isequal(cohortRoot, 0)
            error('extended_script12_sex_stratified_run:NoRoot', 'No cohort folder selected.');
        end
    end

    paths = extended_script8_resolve_paths(cohortRoot);
    need = intersect(paths.missing, ["profilesXlsx"]);
    if ~isempty(need)
        error('extended_script12_sex_stratified_run:MissingInputs', ...
            'Missing Script 7 profiles: %s. Run Script 7 first.', strjoin(need, ', '));
    end

    pal = extended_tol_bright_palette();
    theme = script12_theme_(cfg, pal);
    [outDirs, modeLabel] = script12_make_output_dirs_(cohortRoot, cfg.plotMode);

    logPath = fullfile(outDirs.logs, sprintf('Script12_SexStratified_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() script12_fclose_(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 12: Sex-stratified activity + pre/post-LD amplitude ===\n');
    fprintf('Cohort:  %s\n', paths.cohortRoot);
    fprintf('Output:  %s\n', outDirs.root);
    fprintf('Mode:    %s | dpi=%g | ext=%s\n', cfg.plotMode, theme.dpi, theme.ext);
    script12_log_(LOG, '%s v%s started', SCRIPT_NAME, SCRIPT_VERSION);
    script12_log_(LOG, 'Cohort: %s', paths.cohortRoot);
    script12_log_(LOG, 'profilesXlsx: %s', paths.profilesXlsx);
    script12_log_(LOG, 'Metrics: PreLD + PostLD first-peak − baseline; F vs M only; BH-FDR within ClusterID×Metric');

    data = script12_load_data_(paths, cfg, LOG);
    panels = cfg.script12.panels;
    facets = cfg.script12.facets;
    resolved = script12_resolve_panels_(data.clusterSummary, panels, LOG);

    ampTable = script12_compute_amplitudes_(data, resolved, facets, cfg, LOG);
    if ~isempty(ampTable)
        writetable(ampTable, fullfile(outDirs.tables, 'Sex_Amplitude_PerMouse.csv'));
        script12_log_(LOG, 'Wrote Sex_Amplitude_PerMouse.csv (%d rows)', height(ampTable));
    end

    statsTable = script12_sex_amplitude_stats_(ampTable, facets, cfg, LOG);
    if ~isempty(statsTable)
        writetable(statsTable, fullfile(outDirs.tables, 'Sex_Amplitude_Stats_BH_FDR.csv'));
        script12_log_(LOG, 'Wrote Sex_Amplitude_Stats_BH_FDR.csv (%d rows; Significant_BH=%d)', ...
            height(statsTable), sum(logical(statsTable.Significant_BH)));
    end

    plotFiles = table(string.empty(0,1), string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
        'VariableNames', {'ClusterID','PlotType','SexLayout','FilePath'});
    nOk = 0;

    for i = 1:numel(resolved)
        R = resolved(i);
        face = script12_cluster_face_(R);
        short = sprintf('%s_C%02d', char(regexprep(char(R.BandName), '[^A-Za-z0-9]+', '_')), R.ClusterRank);

        try
            figA = script12_make_activity_grid_(data.activityZT, R, facets, theme, pal, face);
            stemA = sprintf('Sex_Activity_%s', short);
            outA = script12_export_fig_(figA, fullfile(outDirs.figures, stemA), theme);
            close(figA);
            if strlength(outA) > 0
                nOk = nOk + 1;
                plotFiles = [plotFiles; {string(R.ClusterID), "Activity", "Female|Male side-by-side", string(outA)}]; %#ok<AGROW>
                script12_log_(LOG, 'Activity %s: %s', short, outA);
            end
        catch ME
            warning('Script12:Activity', 'Activity %s failed: %s', short, ME.message);
            script12_log_(LOG, 'Activity %s FAILED: %s', short, ME.message);
        end

        ampSpecs = [ ...
            struct('Metric', "PostLD_FirstPeak_au", 'Stem', "PostLD", 'PlotType', "Amplitude_PostLD_FirstPeak"); ...
            struct('Metric', "PreLD_FirstPeak_au", 'Stem', "PreLD", 'PlotType', "Amplitude_PreLD_FirstPeak")];
        for ai = 1:numel(ampSpecs)
            spec = ampSpecs(ai);
            try
                figL = script12_make_amplitude_fig_(ampTable, statsTable, R, facets, theme, pal, face, spec.Metric);
                stemL = sprintf('Sex_Amplitude_%s_%s', spec.Stem, short);
                outL = script12_export_fig_(figL, fullfile(outDirs.figures, stemL), theme);
                close(figL);
                if strlength(outL) > 0
                    nOk = nOk + 1;
                    plotFiles = [plotFiles; {string(R.ClusterID), string(spec.PlotType), ...
                        "F/M overlaid + BH-FDR stars", string(outL)}]; %#ok<AGROW>
                    script12_log_(LOG, 'Amplitude %s %s: %s', spec.Stem, short, outL);
                end
            catch ME
                warning('Script12:Amp', 'Amplitude %s %s failed: %s', spec.Stem, short, ME.message);
                script12_log_(LOG, 'Amplitude %s %s FAILED: %s', spec.Stem, short, ME.message);
            end
        end
    end

    settings = script12_settings_table_(SCRIPT_NAME, SCRIPT_VERSION, paths, cfg, facets);
    xlsxOut = fullfile(outDirs.tables, 'SexStratifiedProfiles_Output.xlsx');
    script12_safe_delete_(xlsxOut);
    script12_safe_writetable_(settings, xlsxOut, 'Settings');
    script12_safe_writetable_(plotFiles, xlsxOut, 'PlotFiles');
    if ~isempty(data.sexSummary)
        script12_safe_writetable_(data.sexSummary, xlsxOut, 'Sex_N_ByCluster');
        writetable(data.sexSummary, fullfile(outDirs.tables, 'Sex_N_ByCluster.csv'));
    end
    if ~isempty(ampTable)
        script12_safe_writetable_(ampTable, xlsxOut, 'Amplitude_PerMouse');
    end
    if ~isempty(statsTable)
        script12_safe_writetable_(statsTable, xlsxOut, 'Amplitude_Stats_BH_FDR');
    end
    writetable(plotFiles, fullfile(outDirs.tables, 'PlotFiles.csv'));

    script12_write_readme_(fullfile(outDirs.root, 'README.md'), resolved, facets, cfg);

    out = struct();
    out.cohortRoot = paths.cohortRoot;
    out.outputRoot = outDirs.root;
    out.modeLabel = modeLabel;
    out.nFigures = nOk;
    out.plotFiles = plotFiles;
    out.amplitudeTable = ampTable;
    out.statsTable = statsTable;
    out.logPath = logPath;

    fprintf('Figures: %d written under %s\n', nOk, outDirs.figures);
    fprintf('Tables:  %s\n', outDirs.tables);
    fprintf('Log:     %s\n', logPath);
    if nOk == 0
        warning('Script12:NoFigures', 'No figures written. Check log: %s', logPath);
    end
end

%% Config / paths
function cfg = script12_fill_cfg_(cfg)
    if ~isfield(cfg, 'script12') || isempty(cfg.script12)
        cfg.script12 = struct();
    end
    if ~isfield(cfg.script12, 'panels') || isempty(cfg.script12.panels)
        cfg.script12.panels = [ ...
            struct('BandName', "UR_1_3", 'ClusterRank', 1); ...
            struct('BandName', "UR_3_6", 'ClusterRank', 1); ...
            struct('BandName', "UR_3_6", 'ClusterRank', 2)];
    end
    if ~isfield(cfg.script12, 'facets') || isempty(cfg.script12.facets)
        cfg.script12.facets = [12, 14, 16, 18, 20, 22];
    end
    % Back-compat aliases from v1.1–1.2
    if (~isfield(cfg.script12, 'ampBaselinePreH') || isempty(cfg.script12.ampBaselinePreH)) && ...
            isfield(cfg.script12, 'ldSpikePreH') && ~isempty(cfg.script12.ldSpikePreH)
        cfg.script12.ampBaselinePreH = double(cfg.script12.ldSpikePreH);
    end
    if ~isfield(cfg.script12, 'ampBaselinePreH') || isempty(cfg.script12.ampBaselinePreH)
        cfg.script12.ampBaselinePreH = 1.0;
    end
    if ~isfield(cfg.script12, 'searchWindowPeriodFrac') || isempty(cfg.script12.searchWindowPeriodFrac)
        cfg.script12.searchWindowPeriodFrac = 1.0;
    end
    if ~isfield(cfg.script12, 'searchWindowCapH') || isempty(cfg.script12.searchWindowCapH)
        if isfield(cfg.script12, 'ldSpikePostH') && ~isempty(cfg.script12.ldSpikePostH)
            cfg.script12.searchWindowCapH = max(3.0, double(cfg.script12.ldSpikePostH));
        else
            cfg.script12.searchWindowCapH = 3.0;
        end
    end
    if ~isfield(cfg.script12, 'searchWindowFloorH') || isempty(cfg.script12.searchWindowFloorH)
        cfg.script12.searchWindowFloorH = 1.0;
    end
    if ~isfield(cfg.script12, 'peakProminenceMin') || isempty(cfg.script12.peakProminenceMin)
        cfg.script12.peakProminenceMin = 0.05;
    end
    if ~isfield(cfg.script12, 'alphaFdr') || isempty(cfg.script12.alphaFdr)
        if isfield(cfg, 'stats') && isfield(cfg.stats, 'alphaFdr') && ~isempty(cfg.stats.alphaFdr)
            cfg.script12.alphaFdr = double(cfg.stats.alphaFdr);
        else
            cfg.script12.alphaFdr = 0.05;
        end
    end
    if ~isfield(cfg.script12, 'minNPerSex') || isempty(cfg.script12.minNPerSex)
        cfg.script12.minNPerSex = 3;
    end
    if ~isfield(cfg.script12, 'bootstrapN') || isempty(cfg.script12.bootstrapN)
        cfg.script12.bootstrapN = 2000;
    end
end

function theme = script12_theme_(cfg, pal)
    theme = struct();
    theme.dpi = cfg.plot.saveDpi;
    theme.ext = cfg.plot.figExt;
    theme.fontName = pal.fontName;
    theme.fontSize = 11;
    theme.labelSize = 12;
    theme.titleSize = 13;
    theme.axesLineWidth = pal.axesLineWidth;
    theme.tickDir = pal.tickDir;
    theme.palette = pal;
end

function [outDirs, modeLabel] = script12_make_output_dirs_(cohortRoot, plotMode)
    switch lower(string(plotMode))
        case "publication"
            modeLabel = "Publication";
        otherwise
            modeLabel = "Development";
    end
    outDirs = struct();
    outDirs.root = fullfile(cohortRoot, sprintf('Script12_SexStratifiedProfiles_%s', char(modeLabel)));
    outDirs.figures = fullfile(outDirs.root, 'Figures');
    outDirs.tables = fullfile(outDirs.root, 'Tables');
    outDirs.logs = fullfile(outDirs.root, 'Logs');
    script12_ensure_dir_(outDirs.root);
    script12_ensure_dir_(outDirs.figures);
    script12_ensure_dir_(outDirs.tables);
    script12_ensure_dir_(outDirs.logs);
end

%% Data
function data = script12_load_data_(paths, cfg, LOG)
    data = struct();
    data.clusterSummary = script12_read_sheet_(paths.profilesXlsx, 'ClusterSummary');
    data.activityZT = script12_read_sheet_(paths.profilesXlsx, 'ActivityComponent_24h');

    if isempty(data.clusterSummary)
        error('extended_script12_sex_stratified_run:NoClusters', ...
            'ClusterSummary missing in %s', paths.profilesXlsx);
    end
    if isempty(data.activityZT)
        error('extended_script12_sex_stratified_run:NoActivity', ...
            'ActivityComponent_24h empty in %s', paths.profilesXlsx);
    end

    data.activityZT = script12_attach_sex_(data.activityZT, LOG, 'Activity');
    data.sexSummary = script12_sex_n_summary_(data.activityZT, cfg.script12.panels, data.clusterSummary);
    data.primaryUR = string(cfg.bands.primaryUR);
end

function T = script12_attach_sex_(T, LOG, tag)
    if isempty(T), return; end
    if ~ismember('SignalID', T.Properties.VariableNames)
        script12_log_(LOG, '%s: no SignalID — cannot attach Sex', tag);
        T.Sex = repmat("Unknown", height(T), 1);
        return;
    end
    T.Sex = script12_infer_sex_(string(T.SignalID));
    script12_log_(LOG, '%s sex rows: Female=%d Male=%d Unknown=%d', tag, ...
        sum(T.Sex == "Female"), sum(T.Sex == "Male"), sum(T.Sex == "Unknown"));
end

function sx = script12_infer_sex_(signalID)
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

function resolved = script12_resolve_panels_(CS, panels, LOG)
    resolved = struct('BandName', {}, 'ClusterRank', {}, 'ClusterID', {}, ...
        'PeriodLow_h', {}, 'PeriodHigh_h', {}, 'PeriodCentre_h', {});
    CS.BandName = string(CS.BandName);
    CS.ClusterID = string(CS.ClusterID);
    for i = 1:numel(panels)
        band = string(panels(i).BandName);
        rank = double(panels(i).ClusterRank);
        cid = script12_cluster_by_rank_(CS, band, rank);
        if strlength(cid) == 0
            error('extended_script12_sex_stratified_run:ClusterMissing', ...
                'No cluster for %s rank %d', char(band), rank);
        end
        hit = CS(string(CS.ClusterID) == cid, :);
        script12_log_(LOG, 'Panel %d: %s rank %d -> %s', i, char(band), rank, char(cid));
        resolved(i).BandName = band;
        resolved(i).ClusterRank = rank;
        resolved(i).ClusterID = cid;
        if ismember('PeriodLow_h', hit.Properties.VariableNames)
            resolved(i).PeriodLow_h = double(hit.PeriodLow_h(1));
            resolved(i).PeriodHigh_h = double(hit.PeriodHigh_h(1));
        else
            resolved(i).PeriodLow_h = NaN;
            resolved(i).PeriodHigh_h = NaN;
        end
        if ismember('PeriodCentre_h', hit.Properties.VariableNames) && isfinite(double(hit.PeriodCentre_h(1)))
            resolved(i).PeriodCentre_h = double(hit.PeriodCentre_h(1));
        elseif isfinite(resolved(i).PeriodLow_h) && isfinite(resolved(i).PeriodHigh_h)
            resolved(i).PeriodCentre_h = 0.5 * (resolved(i).PeriodLow_h + resolved(i).PeriodHigh_h);
        else
            resolved(i).PeriodCentre_h = NaN;
        end
    end
end

function cid = script12_cluster_by_rank_(CS, bandName, rank)
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
    end
end

function face = script12_cluster_face_(R)
    if isfinite(R.PeriodLow_h) && isfinite(R.PeriodHigh_h)
        face = sprintf('%s C%02d (%.2f–%.2f h)', char(R.BandName), R.ClusterRank, ...
            R.PeriodLow_h, R.PeriodHigh_h);
    else
        face = sprintf('%s C%02d', char(R.BandName), R.ClusterRank);
    end
end

function S = script12_sex_n_summary_(Act, panels, CS)
    rows = {};
    for i = 1:numel(panels)
        band = string(panels(i).BandName);
        rank = double(panels(i).ClusterRank);
        cid = script12_cluster_by_rank_(CS, band, rank);
        if strlength(cid) == 0, continue; end
        for sx = ["Female", "Male", "Unknown"]
            nAct = script12_n_mice_in_(Act, cid, sx);
            rows(end + 1, :) = {cid, band, rank, sx, nAct}; %#ok<AGROW>
        end
    end
    if isempty(rows)
        S = table();
        return;
    end
    S = cell2table(rows, 'VariableNames', ...
        {'ClusterID','BandName','ClusterRank','Sex','N_Mice_Activity'});
end

function n = script12_n_mice_in_(T, clusterID, sex)
    n = 0;
    if isempty(T) || ~ismember('SignalID', T.Properties.VariableNames), return; end
    B = T;
    if ismember('ClusterID', B.Properties.VariableNames)
        B = B(string(B.ClusterID) == string(clusterID), :);
    end
    if ismember('Sex', B.Properties.VariableNames)
        B = B(string(B.Sex) == string(sex), :);
    end
    if isempty(B), return; end
    n = numel(unique(string(B.SignalID)));
end

%% Amplitude metrics — pre-LD and post-LD first local max
function amp = script12_compute_amplitudes_(data, resolved, facets, cfg, LOG)
    rows = {};
    preH = double(cfg.script12.ampBaselinePreH);
    frac = double(cfg.script12.searchWindowPeriodFrac);
    capH = double(cfg.script12.searchWindowCapH);
    floorH = double(cfg.script12.searchWindowFloorH);
    promMin = double(cfg.script12.peakProminenceMin);
    facets = double(facets(:))';
    metricSpecs = [ ...
        struct('Metric', "PostLD_FirstPeak_au", 'Fn', @script12_first_peak_postld_); ...
        struct('Metric', "PreLD_FirstPeak_au", 'Fn', @script12_first_peak_preld_)];

    Act = data.activityZT;
    neededA = {'ClusterID','Photoperiod_h','SignalID','Sex','ZTBinCenter_h','Activity_zscored'};
    if isempty(Act) || ~all(ismember(neededA, Act.Properties.VariableNames))
        amp = table();
        script12_log_(LOG, 'No amplitude rows: ActivityComponent_24h incomplete');
        return;
    end

    for i = 1:numel(resolved)
        R = resolved(i);
        cid = string(R.ClusterID);
        W = script12_search_window_h_(R, frac, floorH, capH);
        script12_log_(LOG, '%s search window W=%.3g h (period centre=%.3g)', ...
            char(cid), W, double(R.PeriodCentre_h));

        A = Act(string(Act.ClusterID) == cid & ismember(Act.Photoperiod_h, facets), :);
        A = A(A.Sex == "Female" | A.Sex == "Male", :);
        if isempty(A), continue; end

        keys = unique(A(:, {'Photoperiod_h','SignalID','Sex'}), 'rows', 'stable');
        for k = 1:height(keys)
            photo = double(keys.Photoperiod_h(k));
            sig = string(keys.SignalID(k));
            sx = string(keys.Sex(k));
            As = A(A.Photoperiod_h == photo & string(A.SignalID) == sig, :);
            if isempty(As), continue; end
            zt = double(As.ZTBinCenter_h);
            y = double(As.Activity_zscored);
            keep = isfinite(zt) & isfinite(y);
            zt = zt(keep); y = y(keep);
            if isempty(zt), continue; end
            [zt, ord] = sort(zt(:));
            y = y(ord);

            ptp = max(y) - min(y);
            tOff = photo;  % lights-off at ZT = photoperiod hours for L12–L22
            if ~(isfinite(tOff) && tOff > 0 && tOff < 24)
                continue;
            end

            for mi = 1:numel(metricSpecs)
                metric = metricSpecs(mi).Metric;
                peakFn = metricSpecs(mi).Fn;
                ampVal = NaN; peakZT = NaN; baseline = NaN; method = "Skipped";
                [ampVal, peakZT, baseline, method] = peakFn(zt, y, tOff, preH, W, promMin);
                rows(end + 1, :) = {cid, string(R.BandName), double(R.ClusterRank), ...
                    photo, sig, sx, metric, ampVal, peakZT, baseline, W, method, ptp}; %#ok<AGROW>
            end
        end
    end

    if isempty(rows)
        amp = table();
        script12_log_(LOG, 'No amplitude rows computed');
        return;
    end
    amp = cell2table(rows, 'VariableNames', ...
        {'ClusterID','BandName','ClusterRank','Photoperiod_h','SignalID','Sex', ...
         'Metric','Value','PeakZT_h','Baseline_z','SearchWindow_h','PeakMethod','PeakToTrough_z'});
    for mi = 1:numel(metricSpecs)
        m = metricSpecs(mi).Metric;
        sub = amp(string(amp.Metric) == m, :);
        nFirst = sum(sub.PeakMethod == "FirstLocalMax");
        nFall = sum(sub.PeakMethod == "WindowMaxFallback");
        script12_log_(LOG, '%s rows: %d (FirstLocalMax=%d, WindowMaxFallback=%d, other=%d)', ...
            char(m), height(sub), nFirst, nFall, height(sub) - nFirst - nFall);
    end
end

function W = script12_search_window_h_(R, frac, floorH, capH)
    pc = NaN;
    if isfield(R, 'PeriodCentre_h'), pc = double(R.PeriodCentre_h); end
    if ~isfinite(pc) && isfield(R, 'PeriodLow_h') && isfield(R, 'PeriodHigh_h')
        if isfinite(R.PeriodLow_h) && isfinite(R.PeriodHigh_h)
            pc = 0.5 * (double(R.PeriodLow_h) + double(R.PeriodHigh_h));
        end
    end
    if ~isfinite(pc) || pc <= 0, pc = 2.0; end
    W = double(frac) * pc;
    W = max(double(floorH), min(double(capH), W));
end

function [ampVal, peakZT, baseline, method] = script12_first_peak_postld_(zt, y, tOff, preH, W, promMin)
% First local max in [tOff, tOff+W] minus mean in [tOff-preH, tOff). Fallback: window max.
    ampVal = NaN; peakZT = NaN; baseline = NaN; method = "Skipped";
    pre = y(zt >= (tOff - preH) & zt < tOff);
    postMask = zt >= tOff & zt <= (tOff + W);
    ztPost = zt(postMask);
    yPost = y(postMask);
    if isempty(pre) || isempty(yPost)
        return;
    end
    baseline = mean(pre);

    % First local maximum in search window
    peakIdx = [];
    nP = numel(yPost);
    for i = 1:nP
        if nP == 1
            isPeak = true;
        elseif i == 1
            isPeak = yPost(i) > yPost(i + 1);
        elseif i == nP
            isPeak = yPost(i) > yPost(i - 1);
        else
            isPeak = (yPost(i) >= yPost(i - 1)) && (yPost(i) > yPost(i + 1));
        end
        if isPeak && (yPost(i) - baseline) >= promMin
            peakIdx = i;
            break;
        end
    end

    if ~isempty(peakIdx)
        ampVal = yPost(peakIdx) - baseline;
        peakZT = ztPost(peakIdx);
        method = "FirstLocalMax";
        return;
    end

    % Fallback: max in search window − baseline
    [ymax, ix] = max(yPost);
    ampVal = ymax - baseline;
    peakZT = ztPost(ix);
    method = "WindowMaxFallback";
end

function [ampVal, peakZT, baseline, method] = script12_first_peak_preld_(zt, y, tOff, preH, W, promMin)
% First local max in [tOff−W, tOff) minus mean in [tOff−W−pre, tOff−W). Fallback: window max.
    ampVal = NaN; peakZT = NaN; baseline = NaN; method = "Skipped";
    tLo = tOff - W;
    preLo = tLo - preH;
    pre = y(zt >= preLo & zt < tLo);
    searchMask = zt >= tLo & zt < tOff;
    ztPre = zt(searchMask);
    yPre = y(searchMask);
    if isempty(pre) || isempty(yPre)
        return;
    end
    baseline = mean(pre);

    peakIdx = [];
    nP = numel(yPre);
    for i = 1:nP
        if nP == 1
            isPeak = true;
        elseif i == 1
            isPeak = yPre(i) > yPre(i + 1);
        elseif i == nP
            isPeak = yPre(i) > yPre(i - 1);
        else
            isPeak = (yPre(i) >= yPre(i - 1)) && (yPre(i) > yPre(i + 1));
        end
        if isPeak && (yPre(i) - baseline) >= promMin
            peakIdx = i;
            break;
        end
    end

    if ~isempty(peakIdx)
        ampVal = yPre(peakIdx) - baseline;
        peakZT = ztPre(peakIdx);
        method = "FirstLocalMax";
        return;
    end

    [ymax, ix] = max(yPre);
    ampVal = ymax - baseline;
    peakZT = ztPre(ix);
    method = "WindowMaxFallback";
end

%% Activity grid (F|M side-by-side)
function fig = script12_make_activity_grid_(T, R, facets, theme, pal, faceLabel)
    facets = double(facets(:))';
    assert(numel(facets) == 6, 'Script 12 expects 6 photoperiod facets (L12–L22).');

    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 1900 920]);
    t = tiledlayout(fig, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    sexes = ["Female", "Male"];
    sexCols = {pal.female, pal.male};
    yLim = script12_activity_ylim_(T, R, facets);

    for si = 1:2
        sx = sexes(si);
        meanCol = sexCols{si};
        colOffset = (si - 1) * 3;
        for fi = 1:6
            row = ceil(fi / 3);
            colInHalf = mod(fi - 1, 3) + 1;
            ax = nexttile(t, (row - 1) * 6 + colOffset + colInHalf);
            hold(ax, 'on'); set(ax, 'Color', 'w');
            photo = facets(fi);
            has = script12_plot_activity_facet_(ax, T, R, photo, sx, meanCol, theme, yLim);
            title(ax, char(script12_pp_label_(photo)), 'FontWeight', 'bold', ...
                'FontName', theme.fontName, 'FontSize', theme.titleSize, 'Interpreter', 'none');
            if ~has
                text(ax, 0.5, 0.5, sprintf('No %s data', lower(char(sx))), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                    'FontName', theme.fontName, 'FontSize', 9, 'Color', [0.45 0.45 0.45]);
            end
            if colInHalf == 1
                ylabel(ax, 'Activity (z)', 'FontWeight', 'bold', 'FontName', theme.fontName);
            end
            if row == 2
                xlabel(ax, 'ZT (h)', 'FontWeight', 'bold', 'FontName', theme.fontName);
            end
            xlim(ax, [0 24]);
            if ~isempty(yLim), ylim(ax, yLim); end
            script12_style_axes_(ax, theme);
        end
    end

    annotation(fig, 'textbox', [0.08 0.955 0.40 0.035], 'String', 'Female', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 14, 'Color', pal.female);
    annotation(fig, 'textbox', [0.52 0.955 0.40 0.035], 'String', 'Male', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontName', theme.fontName, 'FontWeight', 'bold', 'FontSize', 14, 'Color', pal.male);
    sgtitle(fig, sprintf('24 h activity (z-scored) — %s  |  Female | Male', faceLabel), ...
        'FontWeight', 'bold', 'FontName', theme.fontName, 'FontSize', 15, 'Interpreter', 'none');
end

function yLim = script12_activity_ylim_(T, R, facets)
    yLim = [-3 3];
    if isempty(T), return; end
    B = T(string(T.ClusterID) == string(R.ClusterID) & ismember(T.Photoperiod_h, facets), :);
    if isempty(B) || ~ismember('Activity_zscored', B.Properties.VariableNames), return; end
    y = double(B.Activity_zscored); y = y(isfinite(y));
    if isempty(y), return; end
    m = max(abs(y));
    pad = ceil(m * 1.1 * 2) / 2;
    if isfinite(pad) && pad >= 1, yLim = [-pad pad]; end
end

function has = script12_plot_activity_facet_(ax, T, R, photo, sex, meanCol, theme, yLim)
    has = false;
    needed = {'ClusterID','Photoperiod_h','SignalID','Sex','ZTBinCenter_h','Activity_zscored'};
    if isempty(T) || ~all(ismember(needed, T.Properties.VariableNames)), return; end
    A = T(string(T.ClusterID) == string(R.ClusterID) & T.Photoperiod_h == photo & string(T.Sex) == string(sex), :);
    if isempty(A), return; end

    sigs = unique(string(A.SignalID), 'stable');
    mouseZT = {}; mouseY = {};
    for s = 1:numel(sigs)
        As = sortrows(A(string(A.SignalID) == sigs(s), :), 'ZTBinCenter_h');
        z = double(As.ZTBinCenter_h); y = double(As.Activity_zscored);
        keep = z >= 0 & z <= 24 & isfinite(y);
        if ~any(keep), continue; end
        plot(ax, z(keep), y(keep), '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        mouseZT{end+1,1} = z(keep); %#ok<AGROW>
        mouseY{end+1,1} = y(keep); %#ok<AGROW>
    end
    if isempty(mouseZT), return; end
    has = true;
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
    plot(ax, centers, meanY, '-', 'Color', meanCol, 'LineWidth', 2.4);
    yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
    if ~isempty(yLim)
        script12_shade_ld_(ax, photo, yLim);
    end
end

%% F vs M amplitude stats (Mann–Whitney + BH-FDR + Cliff's δ)
function S = script12_sex_amplitude_stats_(amp, facets, cfg, LOG)
%SCRIPT12_SEX_AMPLITUDE_STATS_ Per photoperiod Female vs Male Mann–Whitney U.
%   Assumptions: independent mice; two-sided nonparametric location shift;
%   no normality required (small / unequal n). Effect size: Cliff's δ with
%   percentile bootstrap 95% CI. FDR: BH within ClusterID × Metric (α = cfg).
    facets = double(facets(:))';
    minN = max(1, round(double(cfg.script12.minNPerSex)));
    alpha = double(cfg.script12.alphaFdr);
    nBoot = max(200, round(double(cfg.script12.bootstrapN)));
    emptyCols = { ...
        'ClusterID','BandName','ClusterRank','Photoperiod_h','Metric','Test', ...
        'N_Female','N_Male','Median_Female','Median_Male','Mean_Female','Mean_Male', ...
        'PValue_raw','CliffsDelta','CliffsDelta_CI_lo','CliffsDelta_CI_hi', ...
        'FDR_Family','SkippedReason'};
    if isempty(amp)
        S = cell2table(cell(0, numel(emptyCols)), 'VariableNames', emptyCols);
        S = extended_bh_fdr(S, 'PValue_raw', alpha);
        S.StarLabel = strings(0, 1);
        script12_log_(LOG, 'Amplitude stats: empty amp table');
        return;
    end

    metrics = unique(string(amp.Metric), 'stable');
    clusters = unique(amp(:, {'ClusterID','BandName','ClusterRank'}), 'rows', 'stable');
    rows = {};
    rng(20260811);

    for ci = 1:height(clusters)
        cid = string(clusters.ClusterID(ci));
        bname = string(clusters.BandName(ci));
        crank = double(clusters.ClusterRank(ci));
        for mi = 1:numel(metrics)
            metric = metrics(mi);
            fam = string(sprintf('%s|%s', char(cid), char(metric)));
            for fi = 1:numel(facets)
                photo = facets(fi);
                B = amp(string(amp.ClusterID) == cid & string(amp.Metric) == metric & ...
                    amp.Photoperiod_h == photo, :);
                vF = double(B.Value(string(B.Sex) == "Female"));
                vM = double(B.Value(string(B.Sex) == "Male"));
                vF = vF(isfinite(vF));
                vM = vM(isfinite(vM));
                nF = numel(vF); nM = numel(vM);
                medF = NaN; medM = NaN; meanF = NaN; meanM = NaN;
                if nF > 0, medF = median(vF); meanF = mean(vF); end
                if nM > 0, medM = median(vM); meanM = mean(vM); end
                p = NaN; d = NaN; dLo = NaN; dHi = NaN; skip = "";
                if nF < minN || nM < minN
                    skip = sprintf('n<%d per sex', minN);
                else
                    try
                        p = ranksum(vF, vM);  % two-sided Mann–Whitney U
                    catch
                        p = NaN;
                        skip = "ranksum_failed";
                    end
                    [d, dLo, dHi] = script12_cliffs_delta_ci_(vF, vM, nBoot);
                end
                rows(end + 1, :) = {cid, bname, crank, photo, metric, ...
                    "MannWhitney_U_two_sided", nF, nM, medF, medM, meanF, meanM, ...
                    p, d, dLo, dHi, fam, skip}; %#ok<AGROW>
            end
        end
    end

    if isempty(rows)
        S = cell2table(cell(0, numel(emptyCols)), 'VariableNames', emptyCols);
        S = extended_bh_fdr(S, 'PValue_raw', alpha);
        S.StarLabel = strings(0, 1);
        return;
    end
    S = cell2table(rows, 'VariableNames', emptyCols);
    S = extended_bh_fdr(S, 'PValue_raw', alpha);
    labs = strings(height(S), 1);
    for i = 1:height(S)
        labs(i) = string(script12_fdr_star_(S.Q_BH(i), S.Significant_BH(i)));
    end
    S.StarLabel = labs;
    script12_log_(LOG, ['Amplitude stats: Mann–Whitney F vs M; minN/sex=%d; BH-FDR α=%.3g ', ...
        'within ClusterID×Metric; Significant_BH=%d/%d'], ...
        minN, alpha, sum(logical(S.Significant_BH)), height(S));
end

function star = script12_fdr_star_(q, sig)
    star = '';
    if ~(logical(sig) && isfinite(q)), return; end
    if q < 0.001
        star = '***';
    elseif q < 0.01
        star = '**';
    else
        star = '*';  % Significant_BH (Q_BH ≤ α)
    end
end

function [d, lo, hi] = script12_cliffs_delta_ci_(x, y, nBoot)
% Cliff's δ = P(X>Y) − P(X<Y); positive ⇒ Female > Male. Percentile bootstrap CI.
    x = x(:); y = y(:);
    d = script12_cliffs_delta_(x, y);
    lo = NaN; hi = NaN;
    nX = numel(x); nY = numel(y);
    if nX < 1 || nY < 1 || nBoot < 1, return; end
    boot = nan(nBoot, 1);
    for b = 1:nBoot
        xb = x(randi(nX, nX, 1));
        yb = y(randi(nY, nY, 1));
        boot(b) = script12_cliffs_delta_(xb, yb);
    end
    boot = boot(isfinite(boot));
    if isempty(boot), return; end
    lo = prctile(boot, 2.5);
    hi = prctile(boot, 97.5);
end

function d = script12_cliffs_delta_(x, y)
    x = x(:); y = y(:);
    if isempty(x) || isempty(y)
        d = NaN; return;
    end
    n = 0; gt = 0; lt = 0;
    for i = 1:numel(x)
        gt = gt + sum(x(i) > y);
        lt = lt + sum(x(i) < y);
        n = n + numel(y);
    end
    d = (gt - lt) / n;
end

%% Amplitude figure — F/M boxplots side-by-side + BH-FDR stars above each photoperiod pair
function fig = script12_make_amplitude_fig_(amp, statsT, R, facets, theme, pal, faceLabel, metric)
% F vs M box-and-whisker per photoperiod (side-by-side); BH-FDR stars above each pair.
    facets = double(facets(:))';
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 920 540]);
    ax = axes(fig); hold(ax, 'on'); set(ax, 'Color', 'w');

    if nargin < 8 || isempty(metric)
        metric = "PostLD_FirstPeak_au";
    end
    metric = string(metric);
    [yLab, ttlPrefix] = script12_amp_figure_labels_(metric);
    ttl = sprintf('%s — %s', ttlPrefix, faceLabel);

    if isempty(amp)
        text(ax, 0.5, 0.5, 'No amplitude data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', theme.fontName);
        title(ax, ttl, 'FontWeight', 'bold', 'Interpreter', 'none');
        return;
    end

    B = amp(string(amp.ClusterID) == string(R.ClusterID) & string(amp.Metric) == metric, :);
    if isempty(B)
        text(ax, 0.5, 0.5, sprintf('No %s rows', char(metric)), 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', theme.fontName, 'Interpreter', 'none');
        title(ax, ttl, 'FontWeight', 'bold', 'Interpreter', 'none');
        return;
    end

    sexes = ["Female", "Male"];
    cols = {pal.female, pal.male};
    dx = 0.20;
    offsets = [-dx, dx];
    xTick = 1:numel(facets);
    topY = nan(numel(facets), 1);
    boxW = 0.28;

    for fi = 1:numel(facets)
        pairTop = nan(1, 2);
        for si = 1:2
            v = double(B.Value(B.Photoperiod_h == facets(fi) & string(B.Sex) == sexes(si)));
            v = v(isfinite(v));
            x = fi + offsets(si);
            yTop = script12_draw_amp_box_(ax, x, v, cols{si}, boxW);
            if isfinite(yTop)
                pairTop(si) = yTop;
            end
        end
        if any(isfinite(pairTop))
            topY(fi) = max(pairTop, [], 'omitnan');
        end
    end

    St = table();
    if ~isempty(statsT)
        St = statsT(string(statsT.ClusterID) == string(R.ClusterID) & string(statsT.Metric) == metric, :);
    end
    yDataMax = max(topY, [], 'omitnan');
    if ~isfinite(yDataMax), yDataMax = 1; end
    pad = 0.08 * max(yDataMax, eps);
    starCeil = yDataMax;
    for fi = 1:numel(facets)
        if isempty(St), continue; end
        row = St(St.Photoperiod_h == facets(fi), :);
        if isempty(row), continue; end
        if ismember('StarLabel', row.Properties.VariableNames)
            lab = char(string(row.StarLabel(1)));
        else
            lab = script12_fdr_star_(row.Q_BH(1), row.Significant_BH(1));
        end
        if isempty(lab), continue; end
        y0 = topY(fi);
        if ~isfinite(y0), y0 = yDataMax; end
        yStar = y0 + pad;
        text(ax, fi, yStar, lab, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontName', theme.fontName, ...
            'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1], ...
            'HandleVisibility', 'off');
        starCeil = max(starCeil, yStar + 1.5 * pad);
    end
    yTop = max([yDataMax + pad, starCeil], [], 'omitnan');
    if ~isfinite(yTop) || yTop <= 0
        yTop = 1;
    end
    ylim(ax, [0, yTop]);

    set(ax, 'XTick', xTick, 'XTickLabel', arrayfun(@(x) char(script12_pp_label_(x)), facets, 'UniformOutput', false));
    xlim(ax, [0.4, numel(facets) + 0.6]);
    xlabel(ax, 'Photoperiod', 'FontWeight', 'bold', 'FontName', theme.fontName);
    ylabel(ax, yLab, 'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');
    title(ax, ttl, 'FontWeight', 'bold', 'FontName', theme.fontName, 'Interpreter', 'none');

    hF = patch(ax, nan, nan, cols{1}, 'FaceAlpha', 0.55, 'EdgeColor', cols{1} * 0.7, 'DisplayName', 'Female');
    hM = patch(ax, nan, nan, cols{2}, 'FaceAlpha', 0.55, 'EdgeColor', cols{2} * 0.7, 'DisplayName', 'Male');
    legend(ax, [hF, hM], {'Female', 'Male'}, 'Location', 'best', 'Box', 'off', 'FontName', theme.fontName);
    script12_style_axes_(ax, theme);
end

function yTop = script12_draw_amp_box_(ax, xCenter, y, faceColor, halfW)
% Tukey box + whiskers; overlay jittered points. Returns max y used (data or whisker).
    yTop = NaN;
    y = y(:);
    y = y(isfinite(y));
    if isempty(y)
        return;
    end
    if numel(y) == 1
        plot(ax, xCenter, y, 'o', 'Color', faceColor, 'MarkerFaceColor', faceColor, ...
            'MarkerSize', 6, 'HandleVisibility', 'off');
        yTop = y;
        return;
    end

    q1 = prctile(y, 25);
    q2 = prctile(y, 50);
    q3 = prctile(y, 75);
    iqr = q3 - q1;
    if ~isfinite(iqr) || iqr <= 0
        iqr = max(max(y) - min(y), eps);
    end
    lo = max(min(y), q1 - 1.5 * iqr);
    hi = min(max(y), q3 + 1.5 * iqr);
    inWhisk = y(y >= lo & y <= hi);
    if isempty(inWhisk)
        wLo = min(y);
        wHi = max(y);
    else
        wLo = min(inWhisk);
        wHi = max(inWhisk);
    end

    patch(ax, xCenter + [-halfW, halfW, halfW, -halfW], [q1, q1, q3, q3], faceColor, ...
        'FaceAlpha', 0.55, 'EdgeColor', faceColor * 0.7, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax, [xCenter - halfW, xCenter + halfW], [q2, q2], '-', 'Color', faceColor * 0.35, ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(ax, [xCenter, xCenter], [wLo, q1], '-', 'Color', faceColor * 0.7, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax, [xCenter, xCenter], [q3, wHi], '-', 'Color', faceColor * 0.7, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax, [xCenter - halfW * 0.55, xCenter + halfW * 0.55], [wLo, wLo], '-', ...
        'Color', faceColor * 0.7, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax, [xCenter - halfW * 0.55, xCenter + halfW * 0.55], [wHi, wHi], '-', ...
        'Color', faceColor * 0.7, 'LineWidth', 1.0, 'HandleVisibility', 'off');

    rng(abs(round(1e3 * xCenter + sum(y))));
    jitter = (rand(numel(y), 1) - 0.5) * halfW * 0.9;
    scatter(ax, xCenter + jitter, y, 22, 'MarkerFaceColor', faceColor, ...
        'MarkerEdgeColor', faceColor * 0.65, 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
    yTop = max([y(:); wHi], [], 'omitnan');
end

function [yLab, ttlPrefix] = script12_amp_figure_labels_(metric)
    switch string(metric)
        case "PreLD_FirstPeak_au"
            yLab = 'Amplitude (pre-light transition; au)';
            ttlPrefix = 'Pre–lights-off first-peak amplitude';
        otherwise
            yLab = 'Amplitude (post-light transition; au)';
            ttlPrefix = 'Post–lights-off first-peak amplitude';
    end
end

function script12_shade_ld_(ax, photoH, yl)
    photoH = double(photoH);
    if ~isfinite(photoH) || photoH >= 24 || photoH <= 0, return; end
    if nargin < 3 || isempty(yl) || numel(yl) < 2, yl = ylim(ax); end
    h = patch(ax, [photoH 24 24 photoH], [yl(1) yl(1) yl(2) yl(2)], [0.88 0.88 0.88], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    try, uistack(h, 'bottom'); catch, end %#ok<CTCH>
    xline(ax, photoH, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.1, 'HandleVisibility', 'off');
end

%% Export / IO
function outPath = script12_export_fig_(fig, destBaseNoExt, theme)
    outPath = '';
    script12_ensure_dir_(fileparts(destBaseNoExt));
    dpi = round(double(theme.dpi));
    if ~isfinite(dpi) || dpi <= 0, dpi = 150; end
    wantExt = lower(char(string(theme.ext)));
    if ~ismember(wantExt, {'.jpg', '.jpeg', '.png'}), wantExt = '.jpeg'; end
    finalPath = [char(destBaseNoExt) wantExt];
    drawnow;
    try
        exportgraphics(fig, finalPath, 'Resolution', dpi);
        if isfile(finalPath), outPath = string(finalPath); return; end
    catch
    end
    tmpDir = fullfile(tempdir, 'BioRhythms_Script12');
    script12_ensure_dir_(tmpDir);
    tmpBase = fullfile(tmpDir, sprintf('s12_%s', datestr(now, 'HHMMSSFFF')));
    try
        print(fig, tmpBase, '-djpeg', sprintf('-r%d', dpi), '-noui');
        produced = [tmpBase '.jpg'];
        if isfile(produced)
            copyfile(produced, finalPath, 'f');
            try, delete(produced); catch, end
            if isfile(finalPath), outPath = string(finalPath); end
        end
    catch
    end
end

function T = script12_read_sheet_(xlsxPath, sheetName)
    T = table();
    if ~isfile(xlsxPath), return; end
    try
        opts = detectImportOptions(xlsxPath, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
        T = readtable(xlsxPath, opts);
    catch
        T = table();
    end
end

function settings = script12_settings_table_(scriptName, scriptVersion, paths, cfg, facets)
    ampMetricsStr = strjoin([ ...
        "PostLD: first local max in [off, off+W] - mean [off-pre, off); ", ...
        "PreLD: first local max in [off-W, off) - mean [off-W-pre, off-W); ", ...
        "W=clamp(frac*period,[floor,cap]); fallback=window max"], "");
    names = ["ScriptName"; "ScriptVersion"; "RunTimestamp"; "CohortRoot"; ...
        "ProfilesXlsx"; "PlotMode"; "Facets"; ...
        "AmpBaseline_PreH"; "SearchWindow_PeriodFrac"; "SearchWindow_CapH"; "SearchWindow_FloorH"; ...
        "PeakProminenceMin"; "Alpha_FDR"; "MinN_PerSex"; "BootstrapN_CliffsDelta"; ...
        "AmpMetrics"; "SexAmpTest"; "AmpPlotStyle"; "Claim"];
    vals = [string(scriptName); string(scriptVersion); ...
        string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')); ...
        string(paths.cohortRoot); string(paths.profilesXlsx); string(cfg.plotMode); ...
        strjoin(string(facets), ", "); ...
        string(cfg.script12.ampBaselinePreH); string(cfg.script12.searchWindowPeriodFrac); ...
        string(cfg.script12.searchWindowCapH); string(cfg.script12.searchWindowFloorH); ...
        string(cfg.script12.peakProminenceMin); ...
        string(cfg.script12.alphaFdr); string(cfg.script12.minNPerSex); ...
        string(cfg.script12.bootstrapN); ...
        string(ampMetricsStr); ...
        "MannWhitney_U_two_sided Female vs Male only; Cliff_delta_bootCI; BH-FDR within ClusterID|Metric"; ...
        "Box-and-whisker F|M side-by-side per photoperiod; BH-FDR stars above pair"; ...
        "Sex-stratified activity + pre/post-LD first-peak amplitude; F vs M stars; no pre-vs-post test"];
    settings = table(names, vals, 'VariableNames', {'Setting', 'Value'});
end

function script12_write_readme_(path, resolved, facets, cfg)
    fid = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0, return; end
    fprintf(fid, '# Script 12 — Sex-stratified activity + pre/post-LD amplitude\n\n');
    fprintf(fid, 'Female | Male activity grids (L12–L22) and **pre– / post–lights-off first-peak amplitude** (F/M overlaid).\n\n');
    fprintf(fid, 'No phase coherence or ridge-power panels. **No pre-vs-post LD test** — Female vs Male only.\n\n');
    fprintf(fid, '## Amplitude metrics\n\n');
    fprintf(fid, 'From Script 7 band-filtered, z-scored activity (`ActivityComponent_24h`), per mouse × photoperiod:\n\n');
    fprintf(fid, 'Search window W = clamp(%.3g × cluster period centre, [%.3g, %.3g] h).\n\n', ...
        double(cfg.script12.searchWindowPeriodFrac), double(cfg.script12.searchWindowFloorH), ...
        double(cfg.script12.searchWindowCapH));
    fprintf(fid, '### Post-LD (reactive)\n\n');
    fprintf(fid, '1. Baseline = mean in [lights-off − %.3g h, lights-off)\n', double(cfg.script12.ampBaselinePreH));
    fprintf(fid, '2. First local max in [lights-off, lights-off + W] with (peak − baseline) ≥ %.3g\n', ...
        double(cfg.script12.peakProminenceMin));
    fprintf(fid, '3. Fallback: window max − baseline\n');
    fprintf(fid, '4. Y-axis: `Amplitude (post-light transition; au)`\n\n');
    fprintf(fid, '### Pre-LD (anticipatory)\n\n');
    fprintf(fid, '1. Baseline = mean in [lights-off − W − %.3g h, lights-off − W)\n', double(cfg.script12.ampBaselinePreH));
    fprintf(fid, '2. First local max in [lights-off − W, lights-off) with (peak − baseline) ≥ %.3g\n', ...
        double(cfg.script12.peakProminenceMin));
    fprintf(fid, '3. Fallback: window max − baseline\n');
    fprintf(fid, '4. Y-axis: `Amplitude (pre-light transition; au)`\n\n');
    fprintf(fid, '## Sex comparison\n\n');
    fprintf(fid, '- **Test:** two-sided Mann–Whitney U Female vs Male at each photoperiod\n');
    fprintf(fid, '- **Gate:** n ≥ %d per sex\n', round(double(cfg.script12.minNPerSex)));
    fprintf(fid, '- **Effect size:** Cliff''s δ with bootstrap 95%% CI (n=%d)\n', ...
        round(double(cfg.script12.bootstrapN)));
    fprintf(fid, '- **FDR:** Benjamini–Hochberg within each `ClusterID|Metric` (6 photoperiods per cluster per Pre/Post; α = %.3g)\n', ...
        double(cfg.script12.alphaFdr));
    fprintf(fid, '- **Plot:** box-and-whisker F|M side-by-side per photoperiod with jittered points; BH-FDR stars above each pair\n');
    fprintf(fid, '- **Stars:** * / ** / *** from Q_BH when Significant_BH\n');
    fprintf(fid, '- **Tables:** `Sex_Amplitude_PerMouse.csv`, `Sex_Amplitude_Stats_BH_FDR.csv` (n and methods live here / README)\n\n');
    fprintf(fid, '## Clusters\n\n');
    for i = 1:numel(resolved)
        fprintf(fid, '- %s (W ≈ %.2f h)\n', script12_cluster_face_(resolved(i)), ...
            script12_search_window_h_(resolved(i), double(cfg.script12.searchWindowPeriodFrac), ...
            double(cfg.script12.searchWindowFloorH), double(cfg.script12.searchWindowCapH)));
    end
    fprintf(fid, '\n## Photoperiods\n\n%s\n', ...
        strjoin(arrayfun(@(x) char(script12_pp_label_(x)), facets, 'UniformOutput', false), ', '));
    fprintf(fid, '\nAnnotate sparse n (especially UR_3_6 C01 males) when interpreting non-significant cells.\n');
    fclose(fid);
end

function lbl = script12_pp_label_(photo)
    photo = double(photo);
    if photo >= 24, lbl = "LL"; else, lbl = string(sprintf('L%g', photo)); end
end

function script12_style_axes_(ax, theme)
    box(ax, 'off');
    set(ax, 'TickDir', theme.tickDir, 'FontName', theme.fontName, ...
        'FontSize', theme.fontSize, 'LineWidth', theme.axesLineWidth);
end

function script12_ensure_dir_(d)
    if ~isfolder(d), mkdir(d); end
end

function script12_safe_writetable_(T, xlsxPath, sheetName)
    if isempty(T), T = table(); end
    sheetName = char(string(sheetName));
    if strlength(string(sheetName)) > 31
        sheetName = char(extractBefore(string(sheetName), 32));
    end
    try
        writetable(T, xlsxPath, 'Sheet', sheetName);
    catch ME
        warning('Script12:WriteSheet', 'Could not write sheet %s: %s', sheetName, ME.message);
    end
end

function script12_safe_delete_(p)
    if isfile(p)
        try, delete(p); catch, end
    end
end

function script12_log_(LOG, fmt, varargin)
    line = sprintf(fmt, varargin{:});
    fprintf('%s\n', line);
    if ~isempty(LOG) && LOG > 0
        fprintf(LOG, '%s\n', line);
    end
end

function script12_fclose_(LOG)
    if ~isempty(LOG) && LOG > 0, fclose(LOG); end
end
