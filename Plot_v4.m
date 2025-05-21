%% ------------------------------------------------------------------------
% Script: Activity Plotting with Condition Groups, Sex Separation, and
%           Light/Dark Shading
%
% PURPOSE:
% This script reads in an Excel file containing mice activity data and extracts
% the light schedule (e.g., "L12") from the filename. It converts time (hr) to days
% and shades the background for light/dark periods.
%
% Users define one or more condition groups. For each group, a name is entered and one 
% or more activity columns (each representing an individual) are selected. In this updated
% version, each individual’s sex (appended in their name as “-M” or “-F”) is recognized and 
% outputs are separated accordingly. For each sex in a condition group, individual activity 
% plots and a sex-specific group-average plot with standard deviation are produced.
%
% All plots are saved as 600 DPI JPEG images in a dedicated folder per group and sex.
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
darkHours = 24 - lightHours;  % Calculate dark hours

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

prompt = {'Enter number of condition groups:'};
dlgtitle = 'Condition Groups';
dims = [1 50];
numCondStr = inputdlg(prompt, dlgtitle, dims);
if isempty(numCondStr)
    error('No condition group number entered. Aborting.');
end
numGroups = str2double(numCondStr{1});
if isnan(numGroups) || numGroups < 1
    error('Invalid number of condition groups entered. Aborting.');
end

% For each group, get a name and let the user select one or more activity columns.
condGroupNames = cell(numGroups,1);
condGroupColIdx = cell(numGroups,1);
globalSelectedIdx = [];  % To compute unified y-axis scale over all groups
for grp = 1:numGroups
    promptMsg = {sprintf('Enter a name for condition group %d:', grp)};
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
    % The indices returned by listdlg refer to activityVars; add 1 to point to the original table columns.
    selectedIdx = grpCondIdx + 1;  
    condGroupColIdx{grp} = selectedIdx;
    
    globalSelectedIdx = [globalSelectedIdx, selectedIdx];
end
globalSelectedIdx = unique(globalSelectedIdx);

%% --- OUTPUT DIRECTORY SELECTION ---
outputDir = uigetdir('', 'Select Output Folder');
if outputDir == 0
    error('Output folder selection canceled.');
end

%% --- CALCULATE GLOBAL MAXIMUM ACTIVITY FOR UNIFIED Y-AXIS SCALE ---
% Use all selected activity columns from across groups.
activityValues = table2array(activityData(:, globalSelectedIdx));
maxActivityValue = max(activityValues, [], 'all', 'omitnan');  % Global maximum ignoring NaNs
maxYScale = maxActivityValue + 10;  % Add padding

%% --- COLOUR PALETTE ---
lightColour = [240, 228, 66] / 255;   % colour-blind friendly yellow
darkColour  = [186, 186, 186] / 255;    % colour-blind friendly grey

%% --- CALCULATE TOTAL DAYS ---
totalDays = ceil(max(timeDays));
lightDuration = lightHours / 24;
darkDuration = darkHours / 24;

%% --- LOOP THROUGH EACH CONDITION GROUP ---
for grp = 1:numGroups
    groupName = condGroupNames{grp};
    colsIdx = condGroupColIdx{grp};  % These are indices in activityData
    % Create a dedicated subfolder for this condition group.
    groupFolder = fullfile(outputDir, ['ConditionGroup_' groupName]);
    if ~exist(groupFolder, 'dir')
        mkdir(groupFolder);
    end
    
    % Initialize matrices to hold signals separately for each sex
    groupSignalsMatrixM = [];  % For males
    groupSignalsMatrixF = [];  % For females
    
    %% --- LOOP THROUGH EACH INDIVIDUAL (COLUMN) IN THE GROUP ---
    for i = 1:length(colsIdx)
        colIdx = colsIdx(i);
        mouseName = allVars{colIdx};  % e.g., "Nr2BC12-M" or "Nr2BC12-F"
        
        % Determine sex from the mouse name (assumes ending with '-M' for male, '-F' for female)
        if endsWith(mouseName, '-M')
            sexStr = 'Male';
        elseif endsWith(mouseName, '-F')
            sexStr = 'Female';
        else
            sexStr = 'Unknown';
            warning('Sex could not be determined for %s. Assigned as Unknown.', mouseName);
        end
        
        % Create/verify subfolder within the condition group folder for the given sex.
        sexFolder = fullfile(groupFolder, sexStr);
        if ~exist(sexFolder, 'dir')
            mkdir(sexFolder);
        end
        
        % Extract activity data for this individual.
        activity = activityData.(mouseName);
        if ~isnumeric(activity)
            activity = str2double(activity);
        end
        activity(~isfinite(activity)) = 0;
        
        % Append activity to the appropriate sex group matrix.
        if strcmp(sexStr, 'Male')
            groupSignalsMatrixM = [groupSignalsMatrixM, activity];
        elseif strcmp(sexStr, 'Female')
            groupSignalsMatrixF = [groupSignalsMatrixF, activity];
        end
        
        %% --- PLOT INDIVIDUAL ACTIVITY WITH LIGHT/DARK SHADING ---
        figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
        hold on;
        
        % --- Shade Light/Dark Background ---
        for d = 0:totalDays
            lightBlockStart = d;
            lightBlockEnd = d + lightDuration;
            darkBlockStart = lightBlockEnd;
            darkBlockEnd = darkBlockStart + darkDuration;
            fill([lightBlockStart, lightBlockEnd, lightBlockEnd, lightBlockStart], [0, 0, maxYScale, maxYScale], ...
                lightColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            fill([darkBlockStart, darkBlockEnd, darkBlockEnd, darkBlockStart], [0, 0, maxYScale, maxYScale], ...
                darkColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % --- Plot Activity Data with Vertical Lines ---
        validIdx = ~isnan(activity);
        x = timeDays(validIdx);
        y = activity(validIdx);
        for k = 1:length(x)
            line([x(k) x(k)], [0 y(k)], 'Color', 'k', 'LineWidth', 2);
        end
        
        % --- Axis Settings ---
        set(gca, 'FontName', 'Times New Roman');
        xlabel('Time (days)');
        ylabel('Activity Intensity');
        title(sprintf('Activity Plot: %s | %s | %s', mouseName, groupName, sexStr));
        xticks(0:0.5:totalDays);
        xlim([0 totalDays]);
        ylim([0, maxYScale]);
        set(gca, 'Box', 'off', 'TickDir', 'out');
        hold off;
        
        % --- Save Individual Plot ---
        safeName = strrep(mouseName, ' ', '_');
        outFileName = fullfile(sexFolder, sprintf('%s_%s_%s_ActivityPlot.jpg', groupName, safeName, sexStr));
        print(gcf, outFileName, '-djpeg', '-r600');
        close(gcf);
    end
    
    %% --- PLOT GROUP AVERAGE ACTIVITY WITH STANDARD DEVIATION (By Sex) ---
    % For Males
    if ~isempty(groupSignalsMatrixM)
        meanActivity = nanmean(groupSignalsMatrixM, 2);
        stdActivity = nanstd(groupSignalsMatrixM, 0, 2);
        
        figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
        hold on;
        % --- Shade Light/Dark Background ---
        for d = 0:totalDays
            lightBlockStart = d;
            lightBlockEnd = d + lightDuration;
            darkBlockStart = lightBlockEnd;
            darkBlockEnd = darkBlockStart + darkDuration;
            fill([lightBlockStart, lightBlockEnd, lightBlockEnd, lightBlockStart], [0, 0, maxYScale, maxYScale], ...
                lightColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            fill([darkBlockStart, darkBlockEnd, darkBlockEnd, darkBlockStart], [0, 0, maxYScale, maxYScale], ...
                darkColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % --- Plot Mean and ± Standard Deviation ---
        plot(timeDays, meanActivity, 'k-', 'LineWidth', 2);
        plot(timeDays, meanActivity + stdActivity, 'k:', 'LineWidth', 1);
        plot(timeDays, meanActivity - stdActivity, 'k:', 'LineWidth', 1);
        
        % --- Axis Settings ---
        set(gca, 'FontName', 'Times New Roman');
        xlabel('Time (days)');
        ylabel('Average Activity Intensity');
        title(sprintf('Average Activity Plot with SD | %s | Male', groupName));
        xticks(0:0.5:totalDays);
        xlim([0 totalDays]);
        ylim([0, maxYScale]);
        set(gca, 'Box', 'off', 'TickDir', 'out');
        hold off;
        
        % --- Save Group Average Plot for Males ---
        maleFolder = fullfile(groupFolder, 'Male');  % Already created in the loop if needed
        outFileNameAvg = fullfile(maleFolder, sprintf('%s_Male_Average_ActivityPlot.jpg', groupName));
        print(gcf, outFileNameAvg, '-djpeg', '-r600');
        close(gcf);
    end
    
    % For Females
    if ~isempty(groupSignalsMatrixF)
        meanActivity = nanmean(groupSignalsMatrixF, 2);
        stdActivity = nanstd(groupSignalsMatrixF, 0, 2);
        
        figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
        hold on;
        % --- Shade Light/Dark Background ---
        for d = 0:totalDays
            lightBlockStart = d;
            lightBlockEnd = d + lightDuration;
            darkBlockStart = lightBlockEnd;
            darkBlockEnd = darkBlockStart + darkDuration;
            fill([lightBlockStart, lightBlockEnd, lightBlockEnd, lightBlockStart], [0, 0, maxYScale, maxYScale], ...
                lightColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            fill([darkBlockStart, darkBlockEnd, darkBlockEnd, darkBlockStart], [0, 0, maxYScale, maxYScale], ...
                darkColour, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % --- Plot Mean and ± Standard Deviation ---
        plot(timeDays, meanActivity, 'k-', 'LineWidth', 2);
        plot(timeDays, meanActivity + stdActivity, 'k:', 'LineWidth', 1);
        plot(timeDays, meanActivity - stdActivity, 'k:', 'LineWidth', 1);
        
        % --- Axis Settings ---
        set(gca, 'FontName', 'Times New Roman');
        xlabel('Time (days)');
        ylabel('Average Activity Intensity');
        title(sprintf('Average Activity Plot with SD | %s | Female', groupName));
        xticks(0:0.5:totalDays);
        xlim([0 totalDays]);
        ylim([0, maxYScale]);
        set(gca, 'Box', 'off', 'TickDir', 'out');
        hold off;
        
        % --- Save Group Average Plot for Females ---
        femaleFolder = fullfile(groupFolder, 'Female');  % Already created in the loop if needed
        outFileNameAvg = fullfile(femaleFolder, sprintf('%s_Female_Average_ActivityPlot.jpg', groupName));
        print(gcf, outFileNameAvg, '-djpeg', '-r600');
        close(gcf);
    end
    
    fprintf('Completed processing for condition group "%s".\n', groupName);
end

disp('Activity plotting completed. All individual and group-average plots have been saved by sex.');
