function outDir = extended_plot_ridge_power_takeaway(resyncXlsxOrStruct, cfg)
%EXTENDED_PLOT_RIDGE_POWER_TAKEAWAY Ridge-power takeaway plots from Extended D.
%
%   outDir = extended_plot_ridge_power_takeaway(resyncXlsxOrStruct)
%   outDir = extended_plot_ridge_power_takeaway(resyncXlsxOrStruct, cfg)
%
%   Accepts either:
%     - path to Ultradian_RidgePhase_Resync_Output.xlsx, or
%     - resync struct returned by extended_ridge_resync_run containing outXLSX.
%
%   Outputs under:
%     {resync output root}/RidgePower_Takeaway3_Figures/
%       Figures/*.jpg
%       Figures/*.tif
%       RidgePower_Takeaway3_Summary.xlsx
%
%   Ported from Kent D.5 / D_RidgePeriodGraphing. This helper keeps the
%   publication-style ridge-power takeaway as part of modular E2.

    if nargin < 2 || isempty(cfg)
        if exist('extended_defaults', 'file') == 2
            cfg = extended_defaults();
        else
            cfg = struct();
        end
    end

    inputFile = resolve_input_file_(resyncXlsxOrStruct);
    if ~isfile(inputFile)
        error('extended_plot_ridge_power_takeaway:MissingInput', ...
            'Input workbook not found: %s', inputFile);
    end

    FIG_DPI = get_cfg_(cfg, {'plot', 'saveDpi'}, 600);
    FONT_NAME = 'Times New Roman';
    FONT_SIZE_AX = 18;
    FONT_SIZE_LABEL = 22;
    FONT_SIZE_TITLE = 24;
    LINE_WIDTH_AX = 1.5;

    TRANSITION_ORDER = {'DL','LD','MidLight','MidDark'};
    BAND_ORDER = {'UR_1_3','UR_3_6','UR_6_9','UR_9_12','UR_12_18'};
    HEATMAP_CLIM = [];
    SIG_MARKER = '*';

    filePath = fileparts(inputFile);
    outDir = fullfile(filePath, 'RidgePower_Takeaway3_Figures');
    figDir = fullfile(outDir, 'Figures');
    ensure_dir_(outDir);
    ensure_dir_(figDir);

    fprintf('\nReading ridge-power results from:\n%s\n', inputFile);

    sheetName = 'RidgePowerStats_BH_FDR';
    opts = detectImportOptions(inputFile, 'Sheet', sheetName, ...
        'VariableNamingRule', 'preserve');
    T = readtable(inputFile, opts);
    if isempty(T)
        error('extended_plot_ridge_power_takeaway:EmptySheet', ...
            'Sheet "%s" is empty.', sheetName);
    end

    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);
    reqCols = {'Photoperiod_h','BandName','TransitionType', ...
        'MeanDifference_PostMinusPre','PValue_raw','Q_BH','Significant_BH'};
    for i = 1:numel(reqCols)
        if ~ismember(reqCols{i}, T.Properties.VariableNames)
            error('extended_plot_ridge_power_takeaway:MissingColumn', ...
                'Required column missing: %s', reqCols{i});
        end
    end

    T.BandName = string(T.BandName);
    T.TransitionType = string(T.TransitionType);
    T.Significant_BH_clean = to_logical_(T.Significant_BH);

    validRows = isfinite(T.Photoperiod_h) & ...
        ismember(T.BandName, string(BAND_ORDER)) & ...
        ismember(T.TransitionType, string(TRANSITION_ORDER)) & ...
        isfinite(T.MeanDifference_PostMinusPre);
    T = T(validRows, :);
    if isempty(T)
        error('extended_plot_ridge_power_takeaway:NoValidRows', ...
            'No valid ridge-power rows found after filtering.');
    end

    photoperiods = sort(unique(T.Photoperiod_h(:))');
    photoLabels = strings(size(photoperiods));
    for i = 1:numel(photoperiods)
        if photoperiods(i) >= 24
            photoLabels(i) = "LL";
        else
            photoLabels(i) = "L" + string(photoperiods(i));
        end
    end

    SummaryCounts = table();
    for t = 1:numel(TRANSITION_ORDER)
        tr = string(TRANSITION_ORDER{t});
        idx = T.TransitionType == tr & T.Significant_BH_clean;
        nInc = sum(T.MeanDifference_PostMinusPre(idx) > 0);
        nDec = sum(T.MeanDifference_PostMinusPre(idx) < 0);
        nSig = sum(idx);
        nAll = sum(T.TransitionType == tr);
        SummaryCounts = [SummaryCounts; table(tr, nAll, nSig, nInc, nDec, ...
            'VariableNames', {'TransitionType','N_Tests','N_Significant_BH', ...
            'N_Significant_Increase','N_Significant_Decrease'})]; %#ok<AGROW>
    end

    SummaryMeanAll = groupsummary(T, 'TransitionType', {'mean','median','std'}, ...
        'MeanDifference_PostMinusPre');
    Tsig = T(T.Significant_BH_clean, :);
    if ~isempty(Tsig)
        SummaryMeanSig = groupsummary(Tsig, 'TransitionType', {'mean','median','std'}, ...
            'MeanDifference_PostMinusPre');
    else
        SummaryMeanSig = table();
    end

    summaryFile = fullfile(outDir, 'RidgePower_Takeaway3_Summary.xlsx');
    writetable(T, summaryFile, 'Sheet', 'RidgePowerStats_BH_FDR_Used');
    writetable(SummaryCounts, summaryFile, 'Sheet', 'Significant_Counts');
    writetable(SummaryMeanAll, summaryFile, 'Sheet', 'Mean_AllRows');
    if ~isempty(SummaryMeanSig)
        writetable(SummaryMeanSig, summaryFile, 'Sheet', 'Mean_SignificantRows');
    end

    set(groot, 'defaultAxesFontName', FONT_NAME);
    set(groot, 'defaultTextFontName', FONT_NAME);
    set(groot, 'defaultAxesFontSize', FONT_SIZE_AX);
    set(groot, 'defaultAxesLineWidth', LINE_WIDTH_AX);
    set(groot, 'defaultAxesTickDir', 'out');
    set(groot, 'defaultAxesBox', 'off');

    for t = 1:numel(TRANSITION_ORDER)
        tr = string(TRANSITION_ORDER{t});
        H = nan(numel(BAND_ORDER), numel(photoperiods));
        Sig = false(numel(BAND_ORDER), numel(photoperiods));
        for b = 1:numel(BAND_ORDER)
            band = string(BAND_ORDER{b});
            for p = 1:numel(photoperiods)
                pp = photoperiods(p);
                idx = T.TransitionType == tr & T.BandName == band & T.Photoperiod_h == pp;
                if any(idx)
                    H(b,p) = mean(T.MeanDifference_PostMinusPre(idx), 'omitnan');
                    Sig(b,p) = any(T.Significant_BH_clean(idx));
                end
            end
        end

        f = figure('Color','w','Position',[100 100 1500 850]);
        imagesc(H);
        axis tight;
        colormap(redblue_colormap_(256));
        if isempty(HEATMAP_CLIM)
            maxAbs = max(abs(H(:)), [], 'omitnan');
            if isempty(maxAbs) || isnan(maxAbs) || maxAbs == 0
                maxAbs = 1;
            end
            clim([-maxAbs maxAbs]);
        else
            clim(HEATMAP_CLIM);
        end
        cb = colorbar;
        ylabel(cb, 'Mean ridge power change (post - pre)', ...
            'FontName', FONT_NAME, 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
        xticks(1:numel(photoperiods));
        xticklabels(photoLabels);
        yticks(1:numel(BAND_ORDER));
        yticklabels(strrep(BAND_ORDER, '_', '\_'));
        xlabel('Photoperiod', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
        ylabel('Ultradian band', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
        title(sprintf('Ridge-power change | %s', tr), ...
            'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');
        set(gca, 'TickDir','out', 'LineWidth', LINE_WIDTH_AX, ...
            'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX);
        hold on;
        for b = 1:numel(BAND_ORDER)
            for p = 1:numel(photoperiods)
                if Sig(b,p)
                    text(p, b, SIG_MARKER, 'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', 'FontSize', 28, ...
                        'FontWeight','bold', 'Color','k');
                end
            end
        end
        hold off;
        safeTr = char(regexprep(tr, '[^\w]', '_'));
        save_figure_pair_(f, figDir, sprintf('RidgePower_Heatmap_%s', safeTr), FIG_DPI, cfg);
        close(f);
    end

    f = figure('Color','w','Position',[100 100 1200 750]);
    x = 1:numel(TRANSITION_ORDER);
    incCounts = zeros(size(x));
    decCounts = zeros(size(x));
    for t = 1:numel(TRANSITION_ORDER)
        tr = string(TRANSITION_ORDER{t});
        row = SummaryCounts.TransitionType == tr;
        incCounts(t) = SummaryCounts.N_Significant_Increase(row);
        decCounts(t) = SummaryCounts.N_Significant_Decrease(row);
    end
    bar(x - 0.18, incCounts, 0.35, 'FaceColor', [0.0 0.45 0.70], 'EdgeColor','k');
    hold on;
    bar(x + 0.18, -decCounts, 0.35, 'FaceColor', [0.80 0.40 0.00], 'EdgeColor','k');
    yline(0, 'k--', 'LineWidth', 1.2);
    xticks(x);
    xticklabels(TRANSITION_ORDER);
    ylabel('Number of significant tests', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
    xlabel('Transition / control anchor', 'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
    title('Direction of FDR-significant ridge-power changes', ...
        'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');
    legend({'Post > Pre', 'Post < Pre'}, 'Location','best', 'Box','off');
    set(gca, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX, ...
        'LineWidth', LINE_WIDTH_AX, 'TickDir','out', 'Box','off');
    save_figure_pair_(f, figDir, 'RidgePower_SignificantDirectionCounts', FIG_DPI, cfg);
    close(f);

    f = figure('Color','w','Position',[100 100 1200 750]);
    meanVals = nan(1, numel(TRANSITION_ORDER));
    semVals = nan(1, numel(TRANSITION_ORDER));
    for t = 1:numel(TRANSITION_ORDER)
        tr = string(TRANSITION_ORDER{t});
        vals = T.MeanDifference_PostMinusPre(T.TransitionType == tr);
        meanVals(t) = mean(vals, 'omitnan');
        semVals(t) = std(vals, 'omitnan') ./ sqrt(sum(~isnan(vals)));
    end
    bar(x, meanVals, 0.65, 'FaceColor', [0.35 0.70 0.90], 'EdgeColor','k');
    hold on;
    errorbar(x, meanVals, semVals, 'k', 'LineStyle','none', 'LineWidth', 1.5);
    yline(0, 'k--', 'LineWidth', 1.2);
    xticks(x);
    xticklabels(TRANSITION_ORDER);
    ylabel('Mean ridge power change (post - pre)', ...
        'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
    xlabel('Transition / control anchor', ...
        'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
    title('Mean ridge-power modulation across all tested bands and photoperiods', ...
        'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');
    set(gca, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX, ...
        'LineWidth', LINE_WIDTH_AX, 'TickDir','out', 'Box','off');
    save_figure_pair_(f, figDir, 'RidgePower_MeanPostMinusPre_ByTransition_AllRows', FIG_DPI, cfg);
    close(f);

    if ~isempty(Tsig)
        f = figure('Color','w','Position',[100 100 1200 750]);
        meanSig = nan(1, numel(TRANSITION_ORDER));
        semSig = nan(1, numel(TRANSITION_ORDER));
        for t = 1:numel(TRANSITION_ORDER)
            tr = string(TRANSITION_ORDER{t});
            vals = Tsig.MeanDifference_PostMinusPre(Tsig.TransitionType == tr);
            meanSig(t) = mean(vals, 'omitnan');
            semSig(t) = std(vals, 'omitnan') ./ sqrt(sum(~isnan(vals)));
        end
        bar(x, meanSig, 0.65, 'FaceColor', [0.90 0.60 0.00], 'EdgeColor','k');
        hold on;
        errorbar(x, meanSig, semSig, 'k', 'LineStyle','none', 'LineWidth', 1.5);
        yline(0, 'k--', 'LineWidth', 1.2);
        xticks(x);
        xticklabels(TRANSITION_ORDER);
        ylabel('Mean FDR-significant ridge power change (post - pre)', ...
            'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
        xlabel('Transition / control anchor', ...
            'FontSize', FONT_SIZE_LABEL, 'FontWeight','bold');
        title('FDR-significant ridge-power modulation', ...
            'FontSize', FONT_SIZE_TITLE, 'FontWeight','bold');
        set(gca, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX, ...
            'LineWidth', LINE_WIDTH_AX, 'TickDir','out', 'Box','off');
        save_figure_pair_(f, figDir, ...
            'RidgePower_MeanPostMinusPre_ByTransition_SignificantOnly', FIG_DPI, cfg);
        close(f);
    end

    fprintf('\nDone.\nFigures written to:\n%s\n', figDir);
end

function inputFile = resolve_input_file_(resyncXlsxOrStruct)
    if isstruct(resyncXlsxOrStruct)
        if isfield(resyncXlsxOrStruct, 'outXLSX')
            inputFile = char(string(resyncXlsxOrStruct.outXLSX));
        else
            error('extended_plot_ridge_power_takeaway:MissingOutXLSX', ...
                'Resync struct must contain outXLSX.');
        end
    else
        inputFile = char(string(resyncXlsxOrStruct));
    end
end

function val = get_cfg_(cfg, pathParts, fallback)
    val = fallback;
    cur = cfg;
    for i = 1:numel(pathParts)
        key = pathParts{i};
        if ~isstruct(cur) || ~isfield(cur, key)
            return;
        end
        cur = cur.(key);
    end
    val = cur;
end

function tf = to_logical_(x)
    if islogical(x)
        tf = x(:);
    elseif isnumeric(x)
        tf = x(:) ~= 0;
    else
        s = lower(strtrim(string(x(:))));
        tf = ismember(s, ["true","1","yes","y","pass","passed"]);
    end
end

function save_figure_pair_(figHandle, outDir, baseName, dpi, cfg)
    ext = '.png';
    saveTiff = false;
    if nargin >= 5 && isstruct(cfg) && isfield(cfg, 'plot')
        if isfield(cfg.plot, 'figExt'), ext = char(cfg.plot.figExt); end
        if isfield(cfg.plot, 'saveTiff'), saveTiff = logical(cfg.plot.saveTiff); end
    end
    outFile = fullfile(outDir, [baseName ext]);
    exportgraphics(figHandle, outFile, 'Resolution', dpi);
    fprintf('Saved: %s\n', outFile);
    if saveTiff
        tifFile = fullfile(outDir, [baseName '.tif']);
        exportgraphics(figHandle, tifFile, 'Resolution', dpi);
    end
end

function cmap = redblue_colormap_(n)
    if nargin < 1
        n = 256;
    end
    n1 = floor(n/2);
    n2 = n - n1;
    blue = [0.0 0.35 0.70];
    white = [1.0 1.0 1.0];
    red = [0.80 0.10 0.10];
    c1 = [linspace(blue(1), white(1), n1)', ...
          linspace(blue(2), white(2), n1)', ...
          linspace(blue(3), white(3), n1)'];
    c2 = [linspace(white(1), red(1), n2)', ...
          linspace(white(2), red(2), n2)', ...
          linspace(white(3), red(3), n2)'];
    cmap = [c1; c2];
end

function ensure_dir_(p)
    if ~exist(p, 'dir')
        mkdir(p);
    end
end
