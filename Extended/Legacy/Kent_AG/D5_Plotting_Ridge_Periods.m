%% Plot_RidgePower_Takeaway3_from_ResyncOutput.m
% -------------------------------------------------------------------------
% Purpose
%   Reads Ultradian_RidgePhase_Resync_Output.xlsx from Code D and generates
%   figures supporting the ridge-power modulation takeaway:
%
%   "Validated ultradian ridge power is strongly modulated across the
%   light-dark cycle, tending to decrease after lights-on and increase after
%   lights-off."
%
% Input
%   Ultradian_RidgePhase_Resync_Output.xlsx
%
% Required sheet
%   RidgePowerStats_BH_FDR
%
% Outputs
%   - JPEG figures, 600 DPI
%   - TIFF figures, 600 DPI
%   - RidgePower_Takeaway3_Summary.xlsx
%
% Developed for Isaiah J. Ting photoperiod pipeline
% -------------------------------------------------------------------------

clear; clc; close all;

%% --------------------------- USER SETTINGS -----------------------------

FIG_DPI = 600;
FONT_NAME = 'Times New Roman';
FONT_SIZE_AX = 18;
FONT_SIZE_LABEL = 22;
FONT_SIZE_TITLE = 24;
LINE_WIDTH_AX = 1.5;

ALPHA_FDR = 0.05;

TRANSITION_ORDER = {'DL','LD','MidLight','MidDark'};
BAND_ORDER = {'UR_1_3','UR_3_6','UR_6_9','UR_9_12','UR_12_18'};

% Heatmap colour limit. Leave empty [] for dynamic symmetric scaling.
HEATMAP_CLIM = [];

% Significance marker
SIG_MARKER = '*';

%% ---------------------------- SELECT FILE ------------------------------

[fileName, filePath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, ...
    'Select Ultradian_RidgePhase_Resync_Output.xlsx');

if isequal(fileName,0)
    error('No file selected. Script stopped.');
end

inputFile = fullfile(filePath, fileName);

outDir = fullfile(filePath, 'RidgePower_Takeaway3_Figures');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figDir = fullfile(outDir, 'Figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

fprintf('\nReading ridge-power results from:\n%s\n', inputFile);

%% ----------------------------- READ TABLE -------------------------------

sheetName = 'RidgePowerStats_BH_FDR';

try
    opts = detectImportOptions(inputFile, 'Sheet', sheetName, ...
        'VariableNamingRule', 'preserve');
    T = readtable(inputFile, opts);
catch ME
    error('Could not read sheet "%s". Error: %s', sheetName, ME.message);
end

if isempty(T)
    error('Sheet "%s" is empty.', sheetName);
end

% Convert variable names to valid MATLAB names while preserving mapping.
origNames = T.Properties.VariableNames;
validNames = matlab.lang.makeValidName(origNames);
T.Properties.VariableNames = validNames;

%% -------------------------- REQUIRED COLUMNS ----------------------------

reqCols = { ...
    'Photoperiod_h', ...
    'BandName', ...
    'TransitionType', ...
    'MeanDifference_PostMinusPre', ...
    'PValue_raw', ...
    'Q_BH', ...
    'Significant_BH'};

for i = 1:numel(reqCols)
    if ~ismember(reqCols{i}, T.Properties.VariableNames)
        error('Required column missing: %s', reqCols{i});
    end
end

%% -------------------------- CLEAN VARIABLES -----------------------------

T.BandName = string(T.BandName);
T.TransitionType = string(T.TransitionType);

% Convert significance column robustly.
if islogical(T.Significant_BH)
    sigBH = T.Significant_BH;
elseif isnumeric(T.Significant_BH)
    sigBH = T.Significant_BH ~= 0;
else
    sigBH = strcmpi(string(T.Significant_BH), "true") | ...
            strcmpi(string(T.Significant_BH), "1") | ...
            strcmpi(string(T.Significant_BH), "yes");
end
T.Significant_BH_clean = sigBH;

% Keep valid rows only.
validRows = ~isnan(T.Photoperiod_h) & ...
            ismember(T.BandName, string(BAND_ORDER)) & ...
            ismember(T.TransitionType, string(TRANSITION_ORDER)) & ...
            ~isnan(T.MeanDifference_PostMinusPre);

T = T(validRows, :);

if isempty(T)
    error('No valid ridge-power rows found after filtering.');
end

% Sort photoperiods.
photoperiods = unique(T.Photoperiod_h);
photoperiods = sort(photoperiods(:)');

photoLabels = strings(size(photoperiods));
for i = 1:numel(photoperiods)
    if photoperiods(i) >= 24
        photoLabels(i) = "LL";
    else
        photoLabels(i) = "L" + string(photoperiods(i));
    end
end

%% ---------------------------- SUMMARY TABLES ----------------------------

% Significant increases/decreases by transition.
SummaryCounts = table();

for t = 1:numel(TRANSITION_ORDER)
    tr = string(TRANSITION_ORDER{t});
    idx = T.TransitionType == tr & T.Significant_BH_clean;

    nInc = sum(T.MeanDifference_PostMinusPre(idx) > 0);
    nDec = sum(T.MeanDifference_PostMinusPre(idx) < 0);
    nSig = sum(idx);
    nAll = sum(T.TransitionType == tr);

    SummaryCounts = [SummaryCounts; table( ...
        tr, nAll, nSig, nInc, nDec, ...
        'VariableNames', {'TransitionType','N_Tests','N_Significant_BH', ...
                          'N_Significant_Increase','N_Significant_Decrease'})]; %#ok<AGROW>
end

% Mean post-minus-pre by transition using all rows.
SummaryMeanAll = groupsummary(T, 'TransitionType', {'mean','median','std'}, ...
    'MeanDifference_PostMinusPre');

% Mean post-minus-pre by transition using significant rows only.
Tsig = T(T.Significant_BH_clean, :);
if ~isempty(Tsig)
    SummaryMeanSig = groupsummary(Tsig, 'TransitionType', {'mean','median','std'}, ...
        'MeanDifference_PostMinusPre');
else
    SummaryMeanSig = table();
end

%% ----------------------------- SAVE SUMMARY -----------------------------

summaryFile = fullfile(outDir, 'RidgePower_Takeaway3_Summary.xlsx');

writetable(T, summaryFile, 'Sheet', 'RidgePowerStats_BH_FDR_Used');
writetable(SummaryCounts, summaryFile, 'Sheet', 'Significant_Counts');
writetable(SummaryMeanAll, summaryFile, 'Sheet', 'Mean_AllRows');

if ~isempty(SummaryMeanSig)
    writetable(SummaryMeanSig, summaryFile, 'Sheet', 'Mean_SignificantRows');
end

fprintf('Summary workbook written:\n%s\n', summaryFile);

%% ----------------------------- FIGURE STYLE -----------------------------

set(groot, 'defaultAxesFontName', FONT_NAME);
set(groot, 'defaultTextFontName', FONT_NAME);
set(groot, 'defaultAxesFontSize', FONT_SIZE_AX);
set(groot, 'defaultAxesLineWidth', LINE_WIDTH_AX);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

%% -------------------------- FIGURE 1: HEATMAPS --------------------------
% One heatmap per transition:
% rows = bands
% columns = photoperiods
% values = mean post-minus-pre ridge power

for t = 1:numel(TRANSITION_ORDER)

    tr = string(TRANSITION_ORDER{t});
    H = nan(numel(BAND_ORDER), numel(photoperiods));
    Sig = false(numel(BAND_ORDER), numel(photoperiods));

    for b = 1:numel(BAND_ORDER)
        band = string(BAND_ORDER{b});

        for p = 1:numel(photoperiods)
            pp = photoperiods(p);

            idx = T.TransitionType == tr & ...
                  T.BandName == band & ...
                  T.Photoperiod_h == pp;

            if any(idx)
                H(b,p) = mean(T.MeanDifference_PostMinusPre(idx), 'omitnan');
                Sig(b,p) = any(T.Significant_BH_clean(idx));
            end
        end
    end

    f = figure('Color','w','Position',[100 100 1500 850]);

    imagesc(H);
    axis tight;

    colormap(redblue_colormap(256));

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

    % Add significance markers.
    hold on;
    for b = 1:numel(BAND_ORDER)
        for p = 1:numel(photoperiods)
            if Sig(b,p)
                text(p, b, SIG_MARKER, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'FontSize', 28, ...
                    'FontWeight','bold', ...
                    'Color','k');
            end
        end
    end
    hold off;

    safeTr = char(regexprep(tr, '[^\w]', '_'));
    save_figure_pair(f, figDir, sprintf('RidgePower_Heatmap_%s', safeTr), FIG_DPI);
    close(f);
end

%% -------- FIGURE 2: SIGNIFICANT INCREASE/DECREASE COUNTS BY TRANSITION ---

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

save_figure_pair(f, figDir, 'RidgePower_SignificantDirectionCounts', FIG_DPI);
close(f);

%% ------- FIGURE 3: MEAN POST-MINUS-PRE BY TRANSITION, ALL ROWS ----------

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

save_figure_pair(f, figDir, 'RidgePower_MeanPostMinusPre_ByTransition_AllRows', FIG_DPI);
close(f);

%% ---- FIGURE 4: MEAN SIGNIFICANT POST-MINUS-PRE BY TRANSITION -----------

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

    save_figure_pair(f, figDir, 'RidgePower_MeanPostMinusPre_ByTransition_SignificantOnly', FIG_DPI);
    close(f);
end

fprintf('\nDone.\nFigures written to:\n%s\n', figDir);

%% ---------------------------- LOCAL FUNCTIONS ---------------------------

function save_figure_pair(figHandle, outDir, baseName, dpi)
    jpgFile = fullfile(outDir, [baseName '.jpg']);
    tifFile = fullfile(outDir, [baseName '.tif']);

    exportgraphics(figHandle, jpgFile, 'Resolution', dpi);
    exportgraphics(figHandle, tifFile, 'Resolution', dpi);

    fprintf('Saved: %s\n', jpgFile);
end

function cmap = redblue_colormap(n)
    % Simple diverging blue-white-red colormap.
    % Negative values = blue, positive values = red.

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