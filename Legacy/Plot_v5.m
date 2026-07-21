%% ------------------------------------------------------------------------
% Script: Activity Plotting with Condition Groups and Light/Dark Shading
%
% PURPOSE:
% This script reads in an Excel file containing mice activity data and
% extracts the light schedule (e.g., "L12") from the filename. It converts
% time (hr) to days and shades the background for light/dark periods.
%
% Users define one or more condition groups. For each group, a name is
% entered and one or more activity columns (each representing an
% individual) are selected. For each condition group, the script:
%   - generates individual activity plots for every selected column
%   - generates a group-average plot with standard deviation across all
%     selected individuals.
%
% All plots are saved as 600 DPI JPEG images in a dedicated folder per
% condition group.
%
% CHANGELOG:
% - 2025-11-09: Removed global y axis scaling across all selected columns.
%               Each individual plot now uses a per-mouse max activity
%               (max + 10) for y limits and shading. Group average plots
%               use max(mean ± SD) + 10 to set y limits.
% - 2025-11-09: Removed sex detection based on "-M" and "-F" suffixes and
%               all sex-specific folders and plots. The script now treats
%               all selected columns within a condition group as one pool
%               and produces a single group-average plot.
% ------------------------------------------------------------------------

%% --- Clear Environment ---
clear all; close all; clc;

%% --- FILE SELECTION ---
[activityFileName, activityPath] = uigetfile('*.xlsx', 'Select the Activity Data File');
if isequal(activityFileName, 0)
    error('File selection canceled.');
end
activityFilePath = fullfile(activityPath, activityFileName);

%% --- LOAD ACTIVITY DATA (Preserve Headers) ---
optsAct = detectImportOptions(activityFilePath);
optsAct.VariableNamingRule = 'preserve';
activityData = readtable(activityFilePath, optsAct);
allVars = activityData.Properties.VariableNames;

%% --- EXTRACT LIGHT HOURS FROM FILENAME ---
[~, activityName, ~] = fileparts(activityFileName);
L_str = regexp(activityName, 'L(\d+)', 'tokens');
if isempty(L_str)
    error('Could not detect light schedule (e.g., L12) from filename.');
end
lightHours = str2double(L_str{1}{1});
darkHours  = 24 - lightHours;  % Calculate dark hours

%% --- ACTIVITY TIME TO DAYS ---
% Assumes time column is named 'Time (hr)'
if ismember('Time (hr)', allVars)
    timeHours = activityData.('Time (hr)');
else
    error('The expected time column "Time (hr)" was not found.');
end
timeDays = timeHours / 24;

%% --- SELECT ACTIVITY COLUMNS VIA CONDITION GROUPS ---
% Available activity columns (assumed to be all columns except the time column)
activityVars = allVars(2:end);

prompt   = {'Enter number of condition groups:'};
dlgtitle = 'Condition Groups';
dims     = [1 50];
numCondStr = inputdlg(prompt, dlgtitle, dims);
if isempty(numCondStr)
    error('No condition group number entered. Aborting.');
end
numGroups = str2double(numCondStr{1});
if isnan(numGroups) || numGroups < 1
    error('Invalid number of condition groups entered. Aborting.');
end

% For each group, get a name and let the user select one or more activity columns.
condGroupNames  = cell(numGroups,1);
condGroupColIdx = cell(numGroups,1);
for grp = 1:numGroups
    promptMsg   = {sprintf('Enter a name for condition group %d:', grp)};
    grpNameCell = inputdlg(promptMsg, 'Condition Group Name', [1 50]);
    if isempty(grpNameCell)
        error('No name entered for condition group %d. Aborting.', grp);
    end
    condGroupNames{grp} = grpNameCell{1};
    
    [grpCondIdx, tf] = listdlg('PromptString', ...
        sprintf('Select the activity columns for condition group "%s":', condGroupNames{grp}), ...
        'SelectionMode', 'multiple', 'ListString', activityVars);
    if ~tf
        error('No activity columns selected for condition group "%s". Aborting.', condGroupNames{grp});
    end
    
    % The indices returned by listdlg refer to activityVars; add 1 to point
    % to the original table columns.
    selectedIdx = grpCondIdx + 1;
    condGroupColIdx{grp} = selectedIdx;
end

%% --- OUTPUT DIRECTORY SELECTION ---
outputDir = uigetdir('', 'Select Output Folder');
if outputDir == 0
    error('Output folder selection canceled.');
end

%% --- COLOUR PALETTE ---
lightColour = [240, 228, 66] / 255;   % colour blind friendly yellow
darkColour  = [186, 186, 186] / 255;  % colour blind friendly grey

%% --- CALCULATE TOTAL DAYS ---
totalDays     = ceil(max(timeDays));
lightDuration = lightHours / 24;
darkDuration  = darkHours / 24;

%% --- LOOP THROUGH EACH CONDITION GROUP ---
for grp = 1:numGroups
    groupName = condGroupNames{grp};
    colsIdx   = condGroupColIdx{grp};  % indices in activityData
    
    % Create a dedicated subfolder for this condition group.
    groupFolder = fullfile(outputDir, ['ConditionGroup_' groupName]);
    if ~exist(groupFolder, 'dir')
        mkdir(groupFolder);
    end
    
    % Matrix to hold signals for this condition group (all individuals)
    groupSignalsMatrix = [];
    
    %% --- LOOP THROUGH EACH INDIVIDUAL (COLUMN) IN THE GROUP ---
    for i = 1:length(colsIdx)
        colIdx    = colsIdx(i);
        mouseName = allVars{colIdx};
        
        % Extract activity data for this individual.
        activity = activityData.(mouseName);
        if ~isnumeric(activity)
            activity = str2double(activity);
        end
        activity(~isfinite(activity)) = 0;
        
        % Append activity to the group matrix.
        groupSignalsMatrix = [groupSignalsMatrix, activity]; %#ok<AGROW>
        
        % Per mouse y scale
        thisMax = max(activity, [], 'omitnan');
        if isempty(thisMax) || ~isfinite(thisMax)
            thisMax = 0;
        end
        maxYScale = thisMax + 10;
        
        %% --- PLOT INDIVIDUAL ACTIVITY WITH LIGHT/DARK SHADING ---
        figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
        hold on;
        
        % Shade Light/Dark Background
        for d = 0:totalDays
            lightBlockStart = d;
            lightBlockEnd   = d + lightDuration;
            darkBlockStart  = lightBlockEnd;
            darkBlockEnd    = darkBlockStart + darkDuration;
            
            fill([lightBlockStart, lightBlockEnd, lightBlockEnd, lightBlockStart], ...
                 [0, 0, maxYScale, maxYScale], ...
                 lightColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            fill([darkBlockStart, darkBlockEnd, darkBlockEnd, darkBlockStart], ...
                 [0, 0, maxYScale, maxYScale], ...
                 darkColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % Plot Activity Data with Vertical Lines
        validIdx = ~isnan(activity);
        x = timeDays(validIdx);
        y = activity(validIdx);
        for k = 1:length(x)
            line([x(k) x(k)], [0 y(k)], 'Color', 'k', 'LineWidth', 2);
        end
        
        % Axis Settings
        set(gca, 'FontName', 'Times New Roman');
        xlabel('Time (days)');
        ylabel('Activity Intensity');
        title(sprintf('Activity Plot:  %s  |  %s', mouseName, groupName));
        xticks(0:0.5:totalDays);
        xlim([0 totalDays]);
        ylim([0, maxYScale]);
        set(gca, 'Box', 'off', 'TickDir', 'out');
        hold off;
        
        % Save Individual Plot
        safeName    = strrep(mouseName, ' ', '_');
        outFileName = fullfile(groupFolder, ...
            sprintf('%s_%s_ActivityPlot.jpg', groupName, safeName));
        print(gcf, outFileName, '-djpeg', '-r600');
        close(gcf);
    end
    
    %% --- PLOT GROUP AVERAGE ACTIVITY WITH STANDARD DEVIATION (All animals) ---
    if ~isempty(groupSignalsMatrix)
        meanActivity = nanmean(groupSignalsMatrix, 2);
        stdActivity  = nanstd(groupSignalsMatrix, 0, 2);
        
        % y scale based on mean ± SD
        maxYScale = max(meanActivity + stdActivity, [], 'omitnan') + 10;
        
        figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
        hold on;
        
        % Shade Light/Dark Background
        for d = 0:totalDays
            lightBlockStart = d;
            lightBlockEnd   = d + lightDuration;
            darkBlockStart  = lightBlockEnd;
            darkBlockEnd    = darkBlockStart + darkDuration;
            
            fill([lightBlockStart, lightBlockEnd, lightBlockEnd, lightBlockStart], ...
                 [0, 0, maxYScale, maxYScale], ...
                 lightColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            fill([darkBlockStart, darkBlockEnd, darkBlockEnd, darkBlockStart], ...
                 [0, 0, maxYScale, maxYScale], ...
                 darkColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % Plot Mean and ± Standard Deviation
        plot(timeDays, meanActivity, 'k-', 'LineWidth', 2);
        plot(timeDays, meanActivity + stdActivity, 'k:', 'LineWidth', 1);
        plot(timeDays, meanActivity - stdActivity, 'k:', 'LineWidth', 1);
        
        % Axis Settings
        set(gca, 'FontName', 'Times New Roman');
        xlabel('Time (days)');
        ylabel('Average Activity Intensity');
        title(sprintf('Average Activity Plot with SD  |  %s', groupName));
        xticks(0:0.5:totalDays);
        xlim([0 totalDays]);
        ylim([0, maxYScale]);
        set(gca, 'Box', 'off', 'TickDir', 'out');
        hold off;
        
        % Save Group Average Plot
        outFileNameAvg = fullfile(groupFolder, ...
            sprintf('%s_Group_Average_ActivityPlot.jpg', groupName));
        print(gcf, outFileNameAvg, '-djpeg', '-r600');
        close(gcf);
    end
    
    fprintf('Completed processing for condition group "%s".\n', groupName);
end

disp('Activity plotting completed. All individual and group-average plots have been saved.');