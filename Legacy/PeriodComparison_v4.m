%% ------------------------------------------------------------------------
% Script: Extract and Compare Top Peaks from Multi-Sheet Peak Data
%         with:
%           1) Individual Plots (xlim = [0,25] with 2-unit ticks & optional y-axis break)
%           2) Comparison Plots with a Transformed X-Axis (excluding 0.5–11.5)
%           3) Distribution Histograms per Lighting Condition & Group 
%              (histograms of PeakPeriod_hr by peak rank)
%
% DESCRIPTION:
% 1) Individual Plots:
%    For each selected sheet (named with a trailing photoperiod, e.g. "CTRL_12"),
%    rows with "Average" in SignalID are removed and the condition is parsed (the text
%    before the first underscore in SignalID). Only rows whose parsed condition match
%    user-specified conditions are kept. Then, for each individual (unique SignalID) the
%    top N peaks (PeakRank ≤ topN) are plotted as black circles on an axis fixed to [0,25]
%    with x-ticks every 2. An optional y-axis break is simulated using overlaid axes if the
%    global maximum of PeakValue_log10 is high.
%
% 2) Comparison Plots:
%    For each peak rank (1 to topN), a separate figure is produced. Aggregated data
%    (mean and std of PeakPeriod_hr by condition and photoperiod) are computed.
%    To simulate an x-axis break (excluding the range [0.5,11.5]), photoperiod values
%    ≥ 11.5 are shifted by subtracting 11, with a small horizontal offset added per
%    condition so that points do not overlap. The x-axis is extended to [0,14] so that
%    the last point is visible, with custom x-ticks/labels reflecting the original
%    photoperiods.
%
% 3) Distribution Histograms:
%    For each selected sheet (lighting condition), the data are further split by group
%    (the condition parsed from SignalID). For each group, a figure is produced with
%    a tiled layout arranged in 2 columns (rows determined by ceil(topN/2)).
%    Each subplot corresponds to a peak rank (1 to topN) and shows a histogram of
%    PeakPeriod_hr with bin limits [0,26] and bin width of 1.5. All subplots share the
%    same x-axis scale ([0,26] with ticks every 2) and the same y-axis (set to the
%    maximum bin count across the subplots for that group).
%
% All plots: no grid, no box, tick marks "out", and all text uses 'Interpreter','none'.
%
% CHANGELOG:
% - 2025-11-09: Original version handled sex by inferring "_Male"/"_Female"
%               suffixes in SignalID and produced sex-specific folders and
%               plots (Male/Female/Unknown).
% - 2025-11-09: Removed all sex detection and sex-specific folders. All
%               analyses (individual plots, comparison plots, distributions)
%               are now pooled across individuals irrespective of sex.
%               Folder structure simplified accordingly.
%% ------------------------------------------------------------------------

%% --- Clear Environment ---
clearvars; close all; clc;

%% --- Select the Excel File ---
[excelFileName, excelPath] = uigetfile('*.xlsx', 'Select the Excel File with Peak Data');
if isequal(excelFileName, 0)
    error('No file selected. Script aborted.');
end
excelFilePath = fullfile(excelPath, excelFileName);

%% --- Get Sheet Names and Allow Selection ---
[~, sheetNames] = xlsfinfo(excelFilePath);
if isempty(sheetNames)
    error('No sheets found in the selected Excel file.');
end
[selectedSheetIdx, ok] = listdlg('PromptString', 'Select sheets to process:', ...
    'SelectionMode', 'multiple', 'ListString', sheetNames);
if ~ok || isempty(selectedSheetIdx)
    error('No sheets selected. Aborting.');
end
selectedSheets = sheetNames(selectedSheetIdx);

%% --- Prompt User for Conditions (Group Names) ---
prompt   = {'Enter the number of conditions:'};
dlgtitle = 'Conditions';
dims     = [1,50];
numCondStr = inputdlg(prompt, dlgtitle, dims);
if isempty(numCondStr)
    error('No condition number entered. Aborting.');
end
numUserConds = str2double(strtrim(numCondStr{1}));
if isnan(numUserConds) || numUserConds < 1 || mod(numUserConds,1) ~= 0
    error('Invalid number of conditions (must be a positive integer).');
end

userConditionNames = cell(numUserConds, 1);
for k = 1:numUserConds
    promptMsg    = {sprintf('Enter name for condition %d (as appears in SignalID):', k)};
    condNameCell = inputdlg(promptMsg, 'Condition Name', [1,50]);
    if isempty(condNameCell)
        error('No name entered for condition %d. Aborting.', k);
    end
    userConditionNames{k} = condNameCell{1};
end

%% --- Prompt for Output Folder ---
resultsFolder = uigetdir(pwd, 'Select or create the folder for saving results');
if resultsFolder == 0
    error('No output folder selected. Script aborted.');
end

% Create main output folders.
indivFolder = fullfile(resultsFolder, 'IndividualPlots');
if ~exist(indivFolder, 'dir'), mkdir(indivFolder); end
compFolder = fullfile(resultsFolder, 'ComparisonPlots');
if ~exist(compFolder, 'dir'), mkdir(compFolder); end
distFolder = fullfile(resultsFolder, 'DistributionHistograms');
if ~exist(distFolder, 'dir'), mkdir(distFolder); end

%% --- Parameters ---
topN           = 5;        % Maximum peaks per individual to consider
indivXLim      = [0, 25];  % Fixed x-axis for individual plots
markerSize     = 6;        % Marker size for individual peak points
offsetSpacing  = 0.3;      % Horizontal offset per condition (comparison plots)
compMarkerSize = 5;        % Marker size for comparison plot points
breakThresholdY = 6;       % Threshold for applying y-axis break in individual plots

%% --- Pass 1: Global Stats for Individual Plots (for y-axis break) ---
globalMaxY = 0;  % Maximum PeakValue_log10 across all selected sheets / conditions

for s = 1:length(selectedSheets)
    tempData = readtable(excelFilePath, 'Sheet', selectedSheets{s});
    
    % Basic checks
    reqCols = {'SignalID','PeakRank','PeakPeriod_hr','PeakValue_log10'};
    for i = 1:numel(reqCols)
        if ~ismember(reqCols{i}, tempData.Properties.VariableNames)
            error('Column %s not found in sheet %s.', reqCols{i}, selectedSheets{s});
        end
    end
    
    % Remove averages
    tempData = tempData(~contains(tempData.SignalID, 'Average', 'IgnoreCase', true), :);
    
    % Parse condition (text before first underscore)
    tempData.ConditionParsed = cellfun(@(x) x(1:find(x=='_',1)-1), ...
                                       tempData.SignalID, 'UniformOutput', false);
    
    % Keep only user-specified conditions
    keepIdx = ismember(tempData.ConditionParsed, userConditionNames);
    tempData = tempData(keepIdx, :);
    
    % Restrict to topN peaks
    tempData = tempData(tempData.PeakRank <= topN, :);
    
    if ~isempty(tempData)
        globalMaxY = max(globalMaxY, max(tempData.PeakValue_log10));
    end
end

if globalMaxY < 1
    globalMaxY = 1;
end

if globalMaxY > breakThresholdY
    Y_break   = 0.3 * globalMaxY;  % bottom axis limit
    gap       = 0.1 * globalMaxY;  % gap between bottom and top portions
    Y_top_min = Y_break + gap;     % start of top axis
else
    Y_break   = globalMaxY;
    gap       = 0;
    Y_top_min = globalMaxY;
end

%% --- Initialise Aggregated Data Storage for Comparison ---
% Each row: {Condition, LightDuration, PeakRank, MeanPeakPeriod, StdPeakPeriod, Count}
aggData = [];

%% --- Process Each Selected Sheet (main pass) ---
for s = 1:length(selectedSheets)
    sheetName = selectedSheets{s};
    fprintf('Processing sheet: %s\n', sheetName);
    
    % Parse photoperiod from the sheet name (assume the last token is numeric)
    parts         = strsplit(sheetName, '_');
    lightDuration = str2double(parts{end});
    if isnan(lightDuration)
        lightDuration = 0;
    end
    
    data = readtable(excelFilePath, 'Sheet', sheetName);
    
    % Basic checks
    reqCols = {'SignalID','PeakRank','PeakPeriod_hr','PeakValue_log10'};
    for i = 1:numel(reqCols)
        if ~ismember(reqCols{i}, data.Properties.VariableNames)
            error('Column %s not found in sheet %s.', reqCols{i}, sheetName);
        end
    end
    
    % Remove averages
    data = data(~contains(data.SignalID, 'Average','IgnoreCase',true), :);
    
    % Parse condition (text before first underscore)
    data.ConditionParsed = cellfun(@(x) x(1:find(x=='_',1)-1), ...
                                   data.SignalID, 'UniformOutput', false);
    
    % Keep only user-specified conditions
    keepIdx = ismember(data.ConditionParsed, userConditionNames);
    data    = data(keepIdx, :);
    
    %% ---- Individual Plots (all individuals pooled, no sex) ----
    uniqueInds = unique(data.SignalID);
    for i = 1:length(uniqueInds)
        sigID   = uniqueInds{i};
        indData = data(strcmp(data.SignalID, sigID) & data.PeakRank <= topN, :);
        if isempty(indData)
            continue;
        end
        
        indData = sortrows(indData, 'PeakRank', 'ascend');
        
        figInd = figure('Visible','off');
        if gap > 0
            ax_bottom = axes('Position',[0.13, 0.15, 0.775, 0.75]);
            ax_top    = axes('Position',[0.13, 0.75, 0.775, 0.20]);
        else
            ax_bottom = axes('Position',[0.13, 0.15, 0.775, 0.75]);
        end
        
        % Plot points
        if gap > 0
            axes(ax_bottom);
            hold on;
            idx_bottom = indData.PeakValue_log10 <= Y_break;
            idx_top    = indData.PeakValue_log10 >= Y_top_min;
            
            if any(idx_bottom)
                plot(indData.PeakPeriod_hr(idx_bottom), indData.PeakValue_log10(idx_bottom), ...
                     'ko','MarkerSize', markerSize, 'MarkerFaceColor','k','LineStyle','none');
            end
            
            if any(idx_top)
                axes(ax_top);
                hold on;
                plot(indData.PeakPeriod_hr(idx_top), indData.PeakValue_log10(idx_top), ...
                     'ko','MarkerSize', markerSize, 'MarkerFaceColor','k','LineStyle','none');
                hold off;
            end
        else
            axes(ax_bottom);
            hold on;
            plot(indData.PeakPeriod_hr, indData.PeakValue_log10, ...
                 'ko','MarkerSize', markerSize, 'MarkerFaceColor','k','LineStyle','none');
            hold off;
        end
        
        % Format axes
        if gap > 0
            axes(ax_bottom);
            set(gca, 'XLim', indivXLim, 'XTick', 0:2:25, 'YLim', [0, Y_break], ...
                'Box','off','TickDir','out','FontName','Times New Roman');
            xlabel('Period (hr)','Interpreter','none');
            ylabel('Power (a.u.)','Interpreter','none');
            title(sprintf('Individual: %s (Sheet: %s)', sigID, sheetName), 'Interpreter','none');
            
            axes(ax_top);
            set(gca, 'XLim', indivXLim, 'XTick', [], ...
                'YLim', [Y_top_min, globalMaxY*1.05], ...
                'Box','off','TickDir','out','FontName','Times New Roman');
            axes(ax_bottom);
        else
            axes(ax_bottom);
            set(gca, 'XLim', indivXLim, 'XTick', 0:2:25, ...
                'YLim', [0, globalMaxY*1.05], ...
                'Box','off','TickDir','out','FontName','Times New Roman');
            xlabel('Period (hr)','Interpreter','none');
            ylabel('Power (a.u.)','Interpreter','none');
            title(sprintf('Individual: %s (Sheet: %s)', sigID, sheetName), 'Interpreter','none');
        end
        
        safeSigID  = strrep(sigID, ' ', '_');
        safeSheet  = strrep(sheetName, ' ', '_');
        outFileInd = fullfile(indivFolder, sprintf('%s_%s_Peaks.jpg', safeSheet, safeSigID));
        print(figInd, outFileInd, '-djpeg','-r600');
        close(figInd);
    end
    
    %% ---- Aggregation for Comparison & (later) Distribution Histograms ----
    % For each condition, aggregate by PeakRank
    for c = 1:length(userConditionNames)
        condName = userConditionNames{c};
        condData = data(strcmp(data.ConditionParsed, condName), :);
        if isempty(condData)
            continue;
        end
        
        for r = 1:topN
            rData = condData(condData.PeakRank == r, :);
            if isempty(rData)
                continue;
            end
            meanPeak  = mean(rData.PeakPeriod_hr, 'omitnan');
            stdPeak   = std(rData.PeakPeriod_hr, 'omitnan');
            countVal  = height(rData);
            
            aggData = [aggData; {condName, lightDuration, r, meanPeak, stdPeak, countVal}]; %#ok<AGROW>
        end
    end
end  % End sheet loop

%% --- Convert Aggregated Data to Table ---
if isempty(aggData)
    warning('No aggregated data created. Skipping comparison plots.');
    aggTable = table();
else
    aggTable = cell2table(aggData, 'VariableNames', ...
        {'Condition','LightDuration','PeakRank','MeanPeakPeriod','StdPeakPeriod','Count'});
end

%% --- Determine Global Y-limits for Comparison Plots ---
if ~isempty(aggTable)
    allMinus    = aggTable.MeanPeakPeriod - aggTable.StdPeakPeriod;
    allPlus     = aggTable.MeanPeakPeriod + aggTable.StdPeakPeriod;
    globalYminC = floor(min(allMinus));
    globalYmaxC = ceil(max(allPlus));
    if globalYminC > 0
        globalYminC = 0;
    end
else
    globalYminC = 0;
    globalYmaxC = 1;
end

%% --- Comparison Plots for Each Peak Rank (pooled across sex) ---
if ~isempty(aggTable)
    for r = 1:topN
        figComp = figure('Visible','off');
        ax = axes('Position',[0.13, 0.15, 0.775, 0.75]);
        hold on;
        
        set(ax, 'XLim', [0, 14], 'XGrid','off','YGrid','off','Box','off', ...
                'TickDir','out','FontName','Times New Roman', ...
                'YLim', [globalYminC, globalYmaxC]);
        ylabel('Period (hr)','Interpreter','none');
        title(sprintf('Mean Periods - Peak Rank %d', r), 'Interpreter','none');
        
        rData = aggTable(aggTable.PeakRank == r, :);
        if ~isempty(rData)
            uniqueConds = unique(rData.Condition);
            nConds      = numel(uniqueConds);
            
            for c = 1:nConds
                condName = uniqueConds{c};
                cdData   = rData(strcmp(rData.Condition, condName), :);
                if isempty(cdData)
                    continue;
                end
                [~, sortIdx] = sort(cdData.LightDuration);
                cdData = cdData(sortIdx, :);
                
                offset = (c - (nConds+1)/2) * offsetSpacing;
                xVals  = cdData.LightDuration;
                newX   = xVals;
                for j = 1:length(xVals)
                    if xVals(j) >= 11.5
                        newX(j) = xVals(j) - 11;
                    end
                end
                newX = newX + offset;
                
                errorbar(newX, cdData.MeanPeakPeriod, cdData.StdPeakPeriod, ...
                    'LineStyle','none','Marker','o','MarkerSize', compMarkerSize, ...
                    'LineWidth', 1.5);
            end
            
            % Custom x-ticks/labels reflecting original photoperiods
            rightOriginal = (12:2:24)';  
            rightTrans    = rightOriginal - 11;  
            customTicks   = [0; rightTrans];
            customLabels  = cell(size(customTicks));
            customLabels{1} = '0';
            for iTick = 2:length(customTicks)
                customLabels{iTick} = num2str(rightOriginal(iTick-1));
            end
            set(ax, 'XTick', customTicks, 'XTickLabel', customLabels);
            xlabel('Light duration (hr)','Interpreter','none');
            
            % Draw "break" symbol near transformed zone
            xl = get(ax, 'XLim');
            yl = get(ax, 'YLim');
            hold on;
            gapB = 0.05;
            plot([0.45, 0.5 - gapB], [yl(1)+0.05*(yl(2)-yl(1)), yl(1)], ...
                'k', 'LineWidth', 1);
            plot([0.5 + gapB, 0.55], [yl(1), yl(1)+0.05*(yl(2)-yl(1))], ...
                'k', 'LineWidth', 1);
            hold off;
            
            legend(uniqueConds, 'Location','best','Interpreter','none');
            legend boxoff;
        end
        
        outFileNameComp = fullfile(compFolder, ...
            sprintf('Comparison_PeakRank_%d_All.jpg', r));
        print(figComp, outFileNameComp, '-djpeg','-r600');
        close(figComp);
    end
end

%% --- Distribution Histograms for Each Lighting Condition & Group (pooled) ---
for s = 1:length(selectedSheets)
    sheetName = selectedSheets{s};
    dataSheet = readtable(excelFilePath, 'Sheet', sheetName);
    
    % Remove averages
    dataSheet = dataSheet(~contains(dataSheet.SignalID, 'Average','IgnoreCase',true), :);
    
    % Parse condition
    dataSheet.ConditionParsed = cellfun(@(x) x(1:find(x=='_',1)-1), ...
                                        dataSheet.SignalID, 'UniformOutput', false);
    
    % Restrict to user conditions
    groups = intersect(unique(dataSheet.ConditionParsed), userConditionNames);
    
    for g = 1:length(groups)
        groupName = groups{g};
        groupData = dataSheet(strcmp(dataSheet.ConditionParsed, groupName), :);
        if isempty(groupData)
            continue;
        end
        
        figHist = figure('Visible','off');
        nRows   = ceil(topN / 2);
        t       = tiledlayout(nRows, 2, 'TileSpacing','compact','Padding','compact');
        sgtitle(sprintf('Distribution of Periods for %s, Sheet: %s', ...
            groupName, sheetName), 'Interpreter','none');
        
        % Determine max bin count for a shared y-range across peak ranks
        maxCount = 0;
        for r = 1:topN
            rData = groupData(groupData.PeakRank == r, :);
            if ~isempty(rData)
                [counts, ~] = histcounts(rData.PeakPeriod_hr, ...
                    'BinLimits',[0,26], 'BinWidth',1.5);
                maxCount = max(maxCount, max(counts));
            end
        end
        
        for r = 1:topN
            axTile = nexttile;
            rData  = groupData(groupData.PeakRank == r, :);
            vals   = rData.PeakPeriod_hr;
            
            if isempty(vals)
                text(0.5,0.5, sprintf('Peak Rank %d\n(No data)', r), ...
                    'HorizontalAlignment','center', ...
                    'Interpreter','none','FontName','Times New Roman');
                xlim(axTile, [0,26]);
                xticks(axTile, 0:2:26);
            else
                histogram(axTile, vals, ...
                    'BinLimits', [0,26], 'BinWidth', 1.5, ...
                    'FaceColor', [0.2 0.6 0.5], ...
                    'EdgeColor','none', 'FaceAlpha', 1);
                xlim(axTile, [0,26]);
                xticks(axTile, 0:2:26);
                title(axTile, sprintf('Peak Rank %d', r), 'Interpreter','none');
            end
            
            if r == 1
                ylabel(axTile, 'Frequency', 'Interpreter','none');
            else
                ylabel(axTile, '');
            end
            xlabel(axTile, 'Period (hr)', 'Interpreter','none');
            set(axTile, 'Box','off','TickDir','out','FontName','Times New Roman');
            if maxCount > 0
                ylim(axTile, [0, maxCount]);
            end
            set(axTile, 'XTickLabelRotation', 0); 
        end
        
        outFileNameHist = fullfile(distFolder, ...
            sprintf('Distribution_%s_%s_All.jpg', sheetName, groupName));
        print(figHist, outFileNameHist, '-djpeg','-r600');
        close(figHist);
    end
end

disp('Processing complete. Individual plots, comparison plots, and distribution histograms have been saved.');