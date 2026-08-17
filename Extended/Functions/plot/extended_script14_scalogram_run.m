function out = extended_script14_scalogram_run(cohortRoot, cfg)
%EXTENDED_SCRIPT14_SCALOGRAM_RUN Publication stitched scalograms (RAW + HSub residual, F|M).
%
%   Recomputes group-mean stitched CWT scalograms from Handoff/ + HSub TimeSeries.
%   Does not re-run per-mouse WP_TS. Outputs four single-panel JPEG/PNG figures with
%   identical canvas size, photoperiod labels above the heatmap, and fixed clim
%   (RAW 0–14; HSub residual 0–10; F|M shared within each signal class).
%
%   Outputs: {cohortRoot}/Script14_Scalograms_{Publication|Development}/
%     RAW_Stitched_F_{cohort}.*
%     RAW_Stitched_M_{cohort}.*
%     HSub_Residual_Stitched_F_{cohort}.*
%     HSub_Residual_Stitched_M_{cohort}.*
%     Settings_Script14.csv

    SCRIPT_NAME = 'extended_script14_scalogram_run';
    SCRIPT_VERSION = '1.0.1';

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end
    if ~isfield(cfg, 'plotMode') || isempty(cfg.plotMode)
        cfg.plotMode = 'development';
    end
    cfg.plotMode = lower(strtrim(char(string(cfg.plotMode))));
    cfg = extended_apply_plot_cfg(cfg);
    cfg = script14_fill_cfg_(cfg);

    if nargin < 1 || isempty(cohortRoot)
        cohortRoot = uigetdir(pwd, 'Select cohort results folder (e.g. C57_LP)');
        if isequal(cohortRoot, 0)
            error('extended_script14:NoRoot', 'No cohort folder selected.');
        end
    end
    cohortRoot = char(string(cohortRoot));
    [~, cohortTag] = fileparts(cohortRoot);

    handoffDir = fullfile(cohortRoot, 'Handoff');
    if ~isfolder(handoffDir)
        error('extended_script14:NoHandoff', 'Handoff folder not found: %s', handoffDir);
    end

    theme = script14_theme_(cfg);
    [outDirs, modeLabel] = script14_make_output_dirs_(cohortRoot, cfg.plotMode);

    logPath = fullfile(outDirs.logs, sprintf('Script14_Scalograms_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() script14_fclose_(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 14: Publication stitched scalograms ===\n');
    fprintf('Cohort:  %s\n', cohortRoot);
    fprintf('Output:  %s\n', outDirs.root);
    fprintf('Mode:    %s | dpi=%g | ext=%s\n', cfg.plotMode, theme.scalogram.dpi, theme.ext);
    script14_log_(LOG, '%s v%s started', SCRIPT_NAME, SCRIPT_VERSION);
    script14_log_(LOG, 'Handoff: %s', handoffDir);

    coreCfg = core_defaults();
    entries = across_load_handoff_cohort(handoffDir);
    script14_log_(LOG, 'Photoperiod segments loaded: %d', numel(entries));

    sexGroups = cellstr(string(cfg.script14.sexGroups));
    script14_assert_groups_present_(entries, sexGroups);

    armKey = char(string(cfg.script14.hsubArm));
    armLabel = char(string(cfg.script14.hsubArmLabel));

    %% Stitch + CWT per sex (RAW and HSub residual)
    rawData = struct('group', {}, 'stitch', {}, 'wt', {}, 'periods_hours', {});
    hsubData = struct('group', {}, 'stitch', {}, 'wt', {}, 'periods_hours', {});

    for gi = 1:numel(sexGroups)
        grpName = sexGroups{gi};
        script14_log_(LOG, 'Stitching RAW group %s ...', grpName);

        stitchRaw = extended_script14_stitch_raw_group(entries, grpName, coreCfg);
        if isempty(stitchRaw.signal)
            error('extended_script14:EmptyRawStitch', ...
                'No RAW stitched signal for group "%s".', grpName);
        end
        FB = wavelet_make_filterbank(numel(stitchRaw.signal), coreCfg.samplingMinutes, coreCfg);
        [wtRaw, periods_hours, ~] = wavelet_compute_cwt(stitchRaw.signal, FB);

        rawData(gi).group = grpName; %#ok<AGROW>
        rawData(gi).stitch = stitchRaw;
        rawData(gi).wt = wtRaw;
        rawData(gi).periods_hours = periods_hours;

        script14_log_(LOG, 'Stitching HSub residual (%s) group %s ...', armLabel, grpName);
        stitchH = extended_script14_stitch_hsub_residual_group(entries, grpName, armKey, coreCfg);
        if isempty(stitchH.signal)
            error('extended_script14:EmptyHSubStitch', ...
                'No HSub residual stitched signal for group "%s".', grpName);
        end
        FBh = wavelet_make_filterbank(numel(stitchH.signal), coreCfg.samplingMinutes, coreCfg);
        [wtH, ~, ~] = wavelet_compute_cwt(stitchH.signal, FBh);

        hsubData(gi).group = grpName; %#ok<AGROW>
        hsubData(gi).stitch = stitchH;
        hsubData(gi).wt = wtH;
        hsubData(gi).periods_hours = periods_hours;
    end

    climRaw = double(cfg.script14.climRaw);
    climHsub = double(cfg.script14.climHSub);
    obsRaw = script14_observed_max_(rawData);
    obsHsub = script14_observed_max_(hsubData);
    script14_log_(LOG, 'clim RAW [%g, %g] (fixed; observed max=%g)', climRaw(1), climRaw(2), obsRaw);
    script14_log_(LOG, 'clim HSub residual [%g, %g] (fixed; observed max=%g)', climHsub(1), climHsub(2), obsHsub);
    script14_warn_clim_clip_(obsRaw, climRaw, 'RAW');
    script14_warn_clim_clip_(obsHsub, climHsub, 'HSub residual');

    %% Export figures
    figPaths = strings(0, 1);
    ext = theme.ext;
    if ~startsWith(ext, '.'), ext = ['.' ext]; end

    figDir = outDirs.figures;
    extended_period_gate_ensure_dir(figDir);

    for gi = 1:numel(rawData)
        grpName = rawData(gi).group;
        outRaw = fullfile(figDir, sprintf('RAW_Stitched_%s_%s%s', grpName, cohortTag, ext));
        extended_render_publication_scalogram(rawData(gi).wt, rawData(gi).periods_hours, ...
            rawData(gi).stitch, outRaw, theme, coreCfg, climRaw);
        script14_log_(LOG, 'Wrote %s', outRaw);
        figPaths(end + 1, 1) = string(outRaw); %#ok<AGROW>
    end

    for gi = 1:numel(hsubData)
        grpName = hsubData(gi).group;
        outH = fullfile(figDir, sprintf('HSub_Residual_Stitched_%s_%s%s', grpName, cohortTag, ext));
        extended_render_publication_scalogram(hsubData(gi).wt, hsubData(gi).periods_hours, ...
            hsubData(gi).stitch, outH, theme, coreCfg, climHsub);
        script14_log_(LOG, 'Wrote %s', outH);
        figPaths(end + 1, 1) = string(outH); %#ok<AGROW>
    end

    photoList = arrayfun(@(e) e.lightHours, entries);
    settingsPath = fullfile(outDirs.root, 'Settings_Script14.csv');
    script14_write_settings_(settingsPath, SCRIPT_NAME, SCRIPT_VERSION, cohortRoot, cohortTag, ...
        cfg, theme, climRaw, climHsub, armKey, armLabel, photoList, figPaths);

    out = struct();
    out.cohortRoot = cohortRoot;
    out.cohortTag = cohortTag;
    out.outputRoot = outDirs.root;
    out.modeLabel = modeLabel;
    out.figurePaths = figPaths;
    out.climRaw = climRaw;
    out.climHSubResidual = climHsub;
    out.settingsPath = settingsPath;
    out.logPath = logPath;
    out.nPhotoperiods = numel(entries);

    fprintf('Settings:        %s\n', settingsPath);
    for fi = 1:numel(figPaths)
        fprintf('Figure:          %s\n', figPaths(fi));
    end
    fprintf('Log:             %s\n', logPath);
end

function cfg = script14_fill_cfg_(cfg)
    if ~isfield(cfg, 'script14') || isempty(cfg.script14)
        cfg.script14 = struct();
    end
    s = cfg.script14;
    if ~isfield(s, 'sexGroups') || isempty(s.sexGroups)
        s.sexGroups = ["F", "M"];
    end
    if ~isfield(s, 'hsubArm') || isempty(s.hsubArm)
        s.hsubArm = cfg.hsub.residualArms.SEL_P360;
    end
    if ~isfield(s, 'hsubArmLabel') || isempty(s.hsubArmLabel)
        s.hsubArmLabel = 'SEL_P360';
    end
    if ~isfield(s, 'climMode') || isempty(s.climMode)
        s.climMode = 'fixed';
    end
    if ~isfield(s, 'climRaw') || isempty(s.climRaw) || numel(s.climRaw) ~= 2
        s.climRaw = [0, 14];
    end
    if ~isfield(s, 'climHSub') || isempty(s.climHSub) || numel(s.climHSub) ~= 2
        s.climHSub = [0, 10];
    end
    if ~isfield(s, 'figureSizePx') || isempty(s.figureSizePx)
        s.figureSizePx = [1280, 640];
    end
    cfg.script14 = s;
end

function theme = script14_theme_(cfg)
    theme = plot_config(cfg.plotMode);
    theme = plot_theme_ensure_scalogram(theme);
    theme.fontName = 'Arial';
    theme.scalogramColormap = 'jet';
    theme.figureSizePx = double(cfg.script14.figureSizePx(:))';
    if numel(theme.figureSizePx) ~= 2
        theme.figureSizePx = [1280, 640];
    end
    ext = char(string(cfg.plot.figExt));
    if strlength(string(ext)) == 0
        ext = ['.' char(theme.scalogram.format)];
    elseif ~startsWith(ext, '.')
        ext = ['.' ext];
    end
    theme.ext = ext;
end

function [outDirs, modeLabel] = script14_make_output_dirs_(cohortRoot, plotMode)
    switch lower(string(plotMode))
        case "publication"
            modeLabel = "Publication";
        otherwise
            modeLabel = "Development";
    end
    outDirs = struct();
    outDirs.root = fullfile(cohortRoot, sprintf('Script14_Scalograms_%s', char(modeLabel)));
    outDirs.figures = fullfile(outDirs.root, 'Figures');
    outDirs.logs = fullfile(outDirs.root, 'Logs');
    extended_period_gate_ensure_dir(outDirs.root);
    extended_period_gate_ensure_dir(outDirs.figures);
    extended_period_gate_ensure_dir(outDirs.logs);
end

function mx = script14_observed_max_(dataStruct)
    mx = 0;
    for i = 1:numel(dataStruct)
        v = max(abs(dataStruct(i).wt(:)), [], 'omitnan');
        if isfinite(v)
            mx = max(mx, v);
        end
    end
end

function script14_warn_clim_clip_(obsMax, clim, label)
    if ~(isfinite(obsMax) && obsMax > clim(2))
        return;
    end
    warning('extended_script14:ClimClip', ...
        '%s |CWT| max (%.3f) exceeds fixed clim [%.3g, %.3g]; values above the bar are clipped.', ...
        label, obsMax, clim(1), clim(2));
end

function script14_assert_groups_present_(entries, sexGroups)
    allNames = {};
    for i = 1:numel(entries)
        if ~isfield(entries(i), 'groups') || isempty(entries(i).groups)
            continue;
        end
        for g = 1:numel(entries(i).groups)
            allNames{end + 1} = entries(i).groups(g).name; %#ok<AGROW>
        end
    end
    allNames = unique(allNames, 'stable');
    for gi = 1:numel(sexGroups)
        if ~any(strcmp(allNames, sexGroups{gi}))
            error('extended_script14:MissingSexGroup', ...
                'Sex group "%s" not found in Handoff groups: %s', ...
                sexGroups{gi}, strjoin(allNames, ', '));
        end
    end
end

function script14_write_settings_(outPath, scriptName, scriptVersion, cohortRoot, cohortTag, ...
        cfg, theme, climRaw, climHsub, armKey, armLabel, photoList, figPaths)

    names = ["ScriptName"; "ScriptVersion"; "RunTimestamp"; "CohortRoot"; "CohortTag"; ...
        "PlotMode"; "ExportDpi"; "ExportExt"; "FigureWidthPx"; "FigureHeightPx"; ...
        "ClimMode"; "Clim_RAW_lo"; "Clim_RAW_hi"; "Clim_HSubResidual_lo"; "Clim_HSubResidual_hi"; ...
        "HSubArmKey"; "HSubArmLabel"; "SexGroups"; "Photoperiods_h"; "FigureCount"];
    vals = [string(scriptName); string(scriptVersion); ...
        string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')); ...
        string(cohortRoot); string(cohortTag); string(cfg.plotMode); ...
        string(theme.scalogram.dpi); string(theme.ext); ...
        string(theme.figureSizePx(1)); string(theme.figureSizePx(2)); ...
        string(cfg.script14.climMode); ...
        string(climRaw(1)); string(climRaw(2)); ...
        string(climHsub(1)); string(climHsub(2)); ...
        string(armKey); string(armLabel); ...
        strjoin(string(cfg.script14.sexGroups), '|'); ...
        strjoin(string(photoList), '|'); string(numel(figPaths))];

    T = table(names, vals, 'VariableNames', {'Setting', 'Value'});
    writetable(T, outPath);
    if ~isempty(figPaths)
        figTbl = table(figPaths, 'VariableNames', {'FigurePath'});
        manifestPath = fullfile(fileparts(outPath), 'FigureManifest_Script14.csv');
        writetable(figTbl, manifestPath);
    end
end

function script14_log_(LOG, fmt, varargin)
    line = sprintf(fmt, varargin{:});
    fprintf('%s\n', line);
    if ~isempty(LOG) && LOG > 0
        fprintf(LOG, '%s\n', line);
    end
end

function script14_fclose_(LOG)
    if ~isempty(LOG) && LOG > 0
        fclose(LOG);
    end
end
