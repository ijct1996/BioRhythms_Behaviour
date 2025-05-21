%% ------------------------------------------------------------------------
% Script: Extract and Compare Top Peaks from Multi-Sheet Peak Data
%         with:
%           1) Individual Plots (xlim = [0,25] with 2-unit ticks & optional y-axis break)
%           2) Separate Comparison Plots with a Transformed X-Axis (excluding 0.5–11.5)
%           3) Distribution Histograms per Lighting Condition & Group 
%              (histograms of PeakPeriod_hr by peak rank)
%
% DESCRIPTION:
% 1) Individual Plots:
%    For each selected sheet (named with a trailing photoperiod, e.g. "CTRL_12"),
%    rows with "Average" in SignalID are removed and the condition is parsed (the text
%    before the first underscore in SignalID). Only rows whose parsed condition match
%    user‐specified conditions are kept. Then, for each individual (unique SignalID) the
%    top 5 peaks (PeakRank ≤ 5) are plotted as black circles on an axis fixed to [0,25]
%    with x‐ticks every 2. An optional y‐axis break is simulated using overlaid axes if the
%    global maximum of PeakValue_log10 is high.
%
% 2) Comparison Plots:
%    For each peak rank (1 to 5) and for each sex (Male, Female, Unknown), a separate
%    figure is produced. Aggregated data (mean and std of PeakPeriod_hr by condition,
%    photoperiod, and sex) are computed. To simulate an x‐axis break (excluding the range
%    [0.5,11.5]), photoperiod values ≥ 11.5 are shifted by subtracting 11, with a small
%    horizontal offset added per condition so that points do not overlap. The x‐axis is
%    extended to [0,14] so that the last point is visible, with custom x‐ticks/labels reflecting
%    the original photoperiods.
%
% 3) Distribution Histograms:
%    For each selected sheet (lighting condition), the data are further split by group (the
%    condition parsed from SignalID) and by sex. For each (group, sex) combination, a figure 
%    is produced with a tiled layout arranged in 2 columns (rows determined by ceil(topN/2)).
%    Each subplot corresponds to a peak rank (1 to 5) and shows a histogram of PeakPeriod_hr 
%    with bin limits [0,25] and bin width of 2. All subplots share the same x‐axis scale ([0,25]
%    with ticks every 2) and the same y‐axis (set to the maximum bin count across the subplots for 
%    that (group,sex)). An overall title indicates the lighting condition, group, and sex.
%
% All plots: no grid, no box, tick marks "out", and all text uses 'Interpreter','none'.
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
prompt = {'Enter the number of conditions:'};
dlgtitle = 'Conditions';
dims = [1,50];
numCondStr = inputdlg(prompt, dlgtitle, dims);
if isempty(numCondStr)
    error('No condition number entered. Aborting.');
end
numUserConds = str2double(numCondStr{1});
if isnan(numUserConds) || numUserConds < 1
    error('Invalid number of conditions.');
end
userConditionNames = cell(numUserConds, 1);
for k = 1:numUserConds
    promptMsg = {sprintf('Enter name for condition %d (as appears in SignalID):', k)};
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

% Create sex-specific subfolders for individual plots.
maleIndivFolder = fullfile(indivFolder, 'Male');
if ~exist(maleIndivFolder, 'dir'), mkdir(maleIndivFolder); end
femaleIndivFolder = fullfile(indivFolder, 'Female');
if ~exist(femaleIndivFolder, 'dir'), mkdir(femaleIndivFolder); end
unknownIndivFolder = fullfile(indivFolder, 'Unknown');
if ~exist(unknownIndivFolder, 'dir'), mkdir(unknownIndivFolder); end

%% --- Parameters ---
topN = 5;               % Maximum peaks per individual
indivXLim = [0, 25];    % Fixed x-axis for individual plots
markerSize = 6;         % Marker size for individual peak points
offsetSpacing = 0.3;    % Horizontal offset per condition (comparison plots)
compMarkerSize = 5;     % Marker size for comparison plot points
breakThresholdY = 6;    % Threshold for applying y-axis break in individual plots

%% --- Pass 1: Global Stats for Individual Plots (for y-axis break) ---
globalMaxY = 0;  % Maximum of PeakValue_log10 across all selected sheets
for s = 1:length(selectedSheets)
    tempData = readtable(excelFilePath, 'Sheet', selectedSheets{s});
    % --- Add Sex information ---
    tempData.Sex = cell(height(tempData),1);
    for i = 1:height(tempData)
        if endsWith(tempData.SignalID{i}, '_Male')
            tempData.Sex{i} = 'Male';
        elseif endsWith(tempData.SignalID{i}, '_Female')
            tempData.Sex{i} = 'Female';
        else
            tempData.Sex{i} = 'Unknown';
        end
    end
    tempData = tempData(~contains(tempData.SignalID, 'Average', 'IgnoreCase', true), :);
    % --- Use find(x=='_',1) instead of strfind to obtain only the first underscore ---
    tempData.ConditionParsed = cellfun(@(x) x(1:find(x=='_',1)-1), tempData.SignalID, 'UniformOutput', false);
    keepIdx = ismember(tempData.ConditionParsed, userConditionNames);
    tempData = tempData(keepIdx, :);
    tempData = tempData(tempData.PeakRank <= topN, :);
    if ~isempty(tempData)
        globalMaxY = max(globalMaxY, max(tempData.PeakValue_log10));
    end
end
if globalMaxY < 1, globalMaxY = 1; end

if globalMaxY > breakThresholdY
    Y_break = 0.3 * globalMaxY;  % Bottom portion y-limit for individual plots
    gap = 0.1 * globalMaxY;      % Gap between bottom and top portions
    Y_top_min = Y_break + gap;   % Top portion lower limit
else
    Y_break = globalMaxY;
    gap = 0;
    Y_top_min = globalMaxY;
end

%% --- Initialise Aggregated Data Storage for Comparison ---
% Each row: {Condition, LightDuration, Sex, PeakRank, MeanPeakPeriod, StdPeakPeriod, Count}
aggData = [];

%% --- Initialise Data Storage for Distribution Histograms ---
% For each condition, each peak rank, and each sex, accumulate PeakPeriod_hr values.
sexCategories = {'Male','Female','Unknown'};
peakRankDist = cell(numUserConds, topN, numel(sexCategories));
for i = 1:numUserConds
    for j = 1:topN
        for k = 1:length(sexCategories)
            peakRankDist{i,j,k} = [];
        end
    end
end

%% --- Process Each Selected Sheet ---
for s = 1:length(selectedSheets)
    sheetName = selectedSheets{s};
    fprintf('Processing sheet: %s\n', sheetName);
    
    % Parse photoperiod from the sheet name (assume the last token is numeric)
    parts = strsplit(sheetName, '_');
    lightDuration = str2double(parts{end});
    if isnan(lightDuration)
        lightDuration = 0;
    end
    
    data = readtable(excelFilePath, 'Sheet', sheetName);
    % --- Add Sex information ---
    data.Sex = cell(height(data),1);
    for i = 1:height(data)
        if endsWith(data.SignalID{i}, '_Male')
            data.Sex{i} = 'Male';
        elseif endsWith(data.SignalID{i}, '_Female')
            data.Sex{i} = 'Female';
        else
            data.Sex{i} = 'Unknown';
        end
    end
    
    reqCols = {'SignalID','PeakRank','PeakPeriod_hr','PeakValue_log10'};
    for i = 1:length(reqCols)
        if ~ismember(reqCols{i}, data.Properties.VariableNames)
            error('Column %s not found in sheet %s.', reqCols{i}, sheetName);
        end
    end
    
    data = data(~contains(data.SignalID, 'Average','IgnoreCase',true), :);
    data.ConditionParsed = cellfun(@(x) x(1:find(x=='_',1)-1), data.SignalID, 'UniformOutput', false);
    keepIdx = ismember(data.ConditionParsed, userConditionNames);
    data = data(keepIdx, :);
    
    %% ---- Individual Plots ----
    uniqueInds = unique(data.SignalID);
    for i = 1:length(uniqueInds)
        sigID = uniqueInds{i};
        indData = data(strcmp(data.SignalID, sigID) & data.PeakRank <= topN, :);
        if isempty(indData), continue; end
        indData = sortrows(indData, 'PeakRank', 'ascend');
        
        % Determine sex from the first row (assumes same for all rows of that SignalID)
        curSex = indData.Sex{1};
        switch curSex
            case 'Male'
                outFolder = maleIndivFolder;
            case 'Female'
                outFolder = femaleIndivFolder;
            otherwise
                outFolder = unknownIndivFolder;
        end
        
        figInd = figure('Visible','off');
        if gap > 0
            ax_bottom = axes('Position',[0.13, 0.15, 0.775, 0.75]);
            ax_top = axes('Position',[0.13, 0.75, 0.775, 0.20]);
        else
            ax_bottom = axes('Position',[0.13, 0.15, 0.775, 0.75]);
        end
        
        % Plot individual peak points as black circles.
        if gap > 0
            axes(ax_bottom);
            hold on;
            idx_bottom = indData.PeakValue_log10 <= Y_break;
            idx_top = indData.PeakValue_log10 >= Y_top_min;
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
            plot(indData.PeakPeriod_hr, indData.PeakValue_log10, 'ko','MarkerSize', markerSize, ...
                'MarkerFaceColor','k','LineStyle','none');
            hold off;
        end
        
        % Format individual plot axes.
        if gap > 0
            axes(ax_bottom);
            set(gca, 'XLim', indivXLim, 'XTick', 0:2:25, 'YLim', [0, Y_break], ...
                'Box','off','TickDir','out','FontName','Times New Roman');
            xlabel('Period (hr)','Interpreter','none');
            ylabel('Power (a.u.)','Interpreter','none');
            title(sprintf('Individual: %s (Sheet: %s) - %s', sigID, sheetName, curSex), 'Interpreter','none');
            
            axes(ax_top);
            set(gca, 'XLim', indivXLim, 'XTick', [], 'YLim', [Y_top_min, globalMaxY*1.05], ...
                'Box','off','TickDir','out','FontName','Times New Roman');
            axes(ax_bottom);
        else
            axes(ax_bottom);
            set(gca, 'XLim', indivXLim, 'XTick', 0:2:25, 'YLim', [0, globalMaxY*1.05], ...
                'Box','off','TickDir','out','FontName','Times New Roman');
            xlabel('Period (hr)','Interpreter','none');
            ylabel('Power (a.u.)','Interpreter','none');
            title(sprintf('Individual: %s (Sheet: %s) - %s', sigID, sheetName, curSex), 'Interpreter','none');
        end
        
        safeSigID = strrep(sigID, ' ', '_');
        safeSheet = strrep(sheetName, ' ', '_');
        outFileName = fullfile(outFolder, sprintf('%s_%s_Peaks.jpg', safeSheet, safeSigID));
        print(figInd, outFileName, '-djpeg','-r600');
        close(figInd);
    end
    
    %% ---- Aggregation for Comparison & Distribution Histograms ----
    % For each condition and each sex, aggregate data.
    for c = 1:length(userConditionNames)
        condName = userConditionNames{c};
        for sexIdx = 1:length(sexCategories)
            sexCat = sexCategories{sexIdx};
            condData = data(strcmp(data.ConditionParsed, condName) & strcmp(data.Sex, sexCat), :);
            if isempty(condData)
                continue;
            end
            for r = 1:topN
                rData = condData(condData.PeakRank == r, :);
                if isempty(rData)
                    continue;
                end
                meanPeak = mean(rData.PeakPeriod_hr, 'omitnan');
                stdPeak = std(rData.PeakPeriod_hr, 'omitnan');
                countVal = height(rData);
                aggData = [aggData; {condName, lightDuration, sexCat, r, meanPeak, stdPeak, countVal}];
                
                % Accumulate PeakPeriod_hr values for distribution histograms.
                idxCond = find(strcmp(userConditionNames, condName));
                if ~isempty(idxCond)
                    peakRankDist{idxCond, r, sexIdx} = [peakRankDist{idxCond, r, sexIdx}; rData.PeakPeriod_hr];
                end
            end
        end
    end
end  % End sheet loop

%% --- Convert Aggregated Data to Table ---
aggTable = cell2table(aggData, 'VariableNames', {'Condition','LightDuration','Sex','PeakRank',...
                                                  'MeanPeakPeriod','StdPeakPeriod','Count'});

%% --- Determine Global Y-limits for Comparison Plots ---
allMinus = aggTable.MeanPeakPeriod - aggTable.StdPeakPeriod;
allPlus  = aggTable.MeanPeakPeriod + aggTable.StdPeakPeriod;
globalYminC = floor(min(allMinus));
globalYmaxC = ceil(max(allPlus));
if globalYminC > 0, globalYminC = 0; end

%% --- Comparison Plots for Each Peak Rank and Sex (Separate Figures) ---
for r = 1:topN
    for sexIdx = 1:length(sexCategories)
        sexCat = sexCategories{sexIdx};
        figComp = figure('Visible','off');
        ax = axes('Position',[0.13, 0.15, 0.775, 0.75]);
        hold on;
        
        set(ax, 'XLim', [0, 14], 'XGrid','off','YGrid','off','Box','off',...
            'TickDir','out','FontName','Times New Roman', 'YLim', [globalYminC, globalYmaxC]);
        ylabel('Period (hr)','Interpreter','none');
        title(sprintf('Mean Periods - Peak Rank %d (%s)', r, sexCat), 'Interpreter','none');
        
        condIdx = strcmp(aggTable.Sex, sexCat) & (aggTable.PeakRank == r);
        uniqueConds = unique(aggTable.Condition(condIdx));
        nConds = length(uniqueConds);
        rData = aggTable(aggTable.PeakRank == r & strcmp(aggTable.Sex, sexCat), :);
        for c = 1:nConds
            condName = uniqueConds{c};
            cdData = rData(strcmp(rData.Condition, condName), :);
            if isempty(cdData)
                continue;
            end
            [~, sortIdx] = sort(cdData.LightDuration);
            cdData = cdData(sortIdx, :);
            offset = (c - (nConds+1)/2) * offsetSpacing;
            xVals = cdData.LightDuration;
            newX = xVals;
            for j = 1:length(xVals)
                if xVals(j) >= 11.5
                    newX(j) = xVals(j) - 11;
                end
            end
            newX = newX + offset;
            errorbar(newX, cdData.MeanPeakPeriod, cdData.StdPeakPeriod, ...
                'LineStyle','none','Marker','o','MarkerSize', compMarkerSize, 'LineWidth', 1.5);
        end
        
        rightOriginal = (12:2:24)';  
        rightTrans = rightOriginal - 11;  
        customTicks = [0; rightTrans];
        customLabels = cell(size(customTicks));
        customLabels{1} = '0';
        for iTick = 2:length(customTicks)
            customLabels{iTick} = num2str(rightOriginal(iTick-1));
        end
        set(ax, 'XTick', customTicks, 'XTickLabel', customLabels);
        xlabel('Light duration (hr)','Interpreter','none');
        
        xl = get(ax, 'XLim');
        yl = get(ax, 'YLim');
        hold on;
        gap = 0.05;
        plot([0.45, 0.5 - gap], [yl(1)+0.05*(yl(2)-yl(1)), yl(1)], 'k', 'LineWidth', 1);
        plot([0.5 + gap, 0.55], [yl(1), yl(1)+0.05*(yl(2)-yl(1))], 'k', 'LineWidth', 1);
        hold off;
        
        legend(uniqueConds, 'Location','best','Interpreter','none');
        legend boxoff;
        
        outFileNameComp = fullfile(compFolder, sprintf('Comparison_PeakRank_%d_%s.jpg', r, sexCat));
        print(figComp, outFileNameComp, '-djpeg','-r600');
        close(figComp);
    end
end

%% --- Distribution Histograms for Each Lighting Condition, Group & Sex ---
for s = 1:length(selectedSheets)
    sheetName = selectedSheets{s};
    dataSheet = readtable(excelFilePath, 'Sheet', sheetName);
    dataSheet = dataSheet(~contains(dataSheet.SignalID, 'Average','IgnoreCase',true), :);
    dataSheet.ConditionParsed = cellfun(@(x) x(1:find(x=='_',1)-1), dataSheet.SignalID, 'UniformOutput', false);
    % Add Sex information.
    dataSheet.Sex = cell(height(dataSheet),1);
    for i = 1:height(dataSheet)
        if endsWith(dataSheet.SignalID{i}, '_Male')
            dataSheet.Sex{i} = 'Male';
        elseif endsWith(dataSheet.SignalID{i}, '_Female')
            dataSheet.Sex{i} = 'Female';
        else
            dataSheet.Sex{i} = 'Unknown';
        end
    end
    groups = intersect(unique(dataSheet.ConditionParsed), userConditionNames);
    for g = 1:length(groups)
        groupName = groups{g};
        groupData = dataSheet(strcmp(dataSheet.ConditionParsed, groupName), :);
        if isempty(groupData)
            continue;
        end
        for sexIdx = 1:length(sexCategories)
            sexCat = sexCategories{sexIdx};
            sexData = groupData(strcmp(groupData.Sex, sexCat), :);
            if isempty(sexData)
                continue;
            end
            figHist = figure('Visible','off');
            nRows = ceil(topN / 2);
            t = tiledlayout(nRows, 2, 'TileSpacing','compact','Padding','compact');
            sgtitle(sprintf('Distribution of Periods for %s, Sheet: %s (%s)', groupName, sheetName, sexCat), 'Interpreter','none');
            
            maxCount = 0;
            for r = 1:topN
                rData = sexData(sexData.PeakRank == r, :);
                if ~isempty(rData)
                    [counts, ~] = histcounts(rData.PeakPeriod_hr, 'BinLimits',[0,26], 'BinWidth',1.5);
                    maxCount = max(maxCount, max(counts));
                end
            end
            
            for r = 1:topN
                axTile = nexttile;
                rData = sexData(sexData.PeakRank == r, :);
                vals = rData.PeakPeriod_hr;
                if isempty(vals)
                    text(0.5,0.5, sprintf('Peak Rank %d\n(No data)', r), 'HorizontalAlignment','center',...
                        'Interpreter','none','FontName','Times New Roman');
                    axTile.XTick = 0:2:26;
                    axTile.XLim = [0,26];
                else
                    histogram(axTile, vals, 'BinLimits', [0,26], 'BinWidth', 1.5, 'FaceColor', [0.2 0.6 0.5],...
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
                set(gca, 'XTickLabelRotation', 0); 
            end
            
            outFileNameHist = fullfile(distFolder, sprintf('Distribution_%s_%s_%s.jpg', sheetName, groupName, sexCat));
            print(figHist, outFileNameHist, '-djpeg','-r600');
            close(figHist);
        end
    end
end

disp('Processing complete. Individual plots, comparison plots, and distribution histograms have been saved.');
