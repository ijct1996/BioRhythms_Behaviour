%% ------------------------------------------------------------------------
% Script: Wavelet Analysis of Mice Behaviour Data
%         (Scalograms with Condition Groups)

% PURPOSE:
% Reads an Excel file containing mice behavioural data. After selecting:
%   - A time column (in hours; converted for plotting),
%   - A light condition column (for drawing vertical markers), and
%   - A number of condition groups (each with a user-specified name and
%     the selection of activity columns corresponding to individuals in
%     that condition),

% the script computes continuous wavelet transform scalograms for:
%   (a) each individual (activity column) in the condition, and
%   (b) the condition-specific average signal across those individuals.
%
% Each individual's scalogram and the average scalogram are saved as
% 600 DPI JPEG files in a dedicated subfolder for the condition.
%
% NOTE: Adjust conversion factors, filter parameters, or aesthetic
%       settings as needed.

% CHANGELOG:
% - 2025-11-09: Original version with sex detection based on "-M"/"-F"
%               suffixes and sex-specific subfolders.
% - 2025-11-09: Removed sex detection and sex-specific subfolders. All
%               selected individuals within a condition group are treated
%               as one pool; a single group-average scalogram is computed.
% ------------------------------------------------------------------------

%% --- Clear environment ---
clearvars; close all; clc;

%% --- Select the Input Excel File ---
[fileName, filePath] = uigetfile('*.xlsx', 'Select the source Excel file');
if isequal(fileName, 0)
    error('No file selected. Script aborted.');
end
sourceFile = fullfile(filePath, fileName);

%% --- Read the source data (preserving original column headers) ---
dataTable = readtable(sourceFile, 'VariableNamingRule', 'preserve');
varNames  = dataTable.Properties.VariableNames;

%% --- User Selections: Time and Light Condition Columns ---
% Select the Time Column (in hours)
[timeIdx, tf] = listdlg('PromptString', 'Select the Time Column (in hours)', ...
                        'SelectionMode', 'single', 'ListString', varNames);
if ~tf
    error('Time column not selected. Aborting.');
end

% Select the Light Condition Column (to draw vertical markers)
[lightCondIdx, tf] = listdlg('PromptString', 'Select the Light Condition Column', ...
                             'SelectionMode', 'single', 'ListString', varNames);
if ~tf
    error('Light condition column not selected. Aborting.');
end

%% --- Prompt for Condition Groups ---
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

% For each condition group, get name and activity columns
condGroupNames  = cell(numGroups,1);
condGroupColIdx = cell(numGroups,1);

for grp = 1:numGroups
    % Name for the condition group
    promptMsg   = {sprintf('Enter a name for condition group %d:', grp)};
    grpNameCell = inputdlg(promptMsg, 'Condition Group Name', [1 50]);
    if isempty(grpNameCell)
        error('No name entered for condition group %d. Aborting.', grp);
    end
    condGroupNames{grp} = grpNameCell{1};
    
    % Select one or more activity columns for this condition
    [grpCondIdx, tf] = listdlg('PromptString', ...
        sprintf('Select the activity data columns for condition group "%s":', condGroupNames{grp}), ...
        'SelectionMode', 'multiple', 'ListString', varNames);
    if ~tf
        error('No activity data columns selected for condition group "%s". Aborting.', condGroupNames{grp});
    end
    condGroupColIdx{grp} = grpCondIdx;
end

%% --- Prompt User for Output Folder ---
resultsFolder = uigetdir(pwd, 'Select or create the folder for saving results');
if resultsFolder == 0
    error('No output folder selected. Script aborted.');
end

%% --- Process Global Time and Light Condition Vectors ---
time_hours = dataTable{:, timeIdx};
time_min   = time_hours * 60;                  % hours -> minutes
time_day   = time_min / (60 * 24);             % minutes -> days
samplingPeriod = mean(diff(time_min));         % sampling period in minutes

lightConditionVector = dataTable{:, lightCondIdx};
condChangeIdx        = find(diff(lightConditionVector) ~= 0);

%% --- Process Each Condition Group ---
for grp = 1:numGroups
    groupName   = condGroupNames{grp};
    activityCols = condGroupColIdx{grp};
    
    % Subfolder for this condition group
    groupFolder = fullfile(resultsFolder, ['Condition_' groupName]);
    if ~exist(groupFolder, 'dir')
        mkdir(groupFolder);
    end
    
    % Optional: subfolder for individual scalograms
    indivFolder = fullfile(groupFolder, 'Individuals');
    if ~exist(indivFolder, 'dir')
        mkdir(indivFolder);
    end
    
    % Matrix to store signals for group average
    signalsMatrix = [];
    
    %% --- Loop over Individuals in this Condition ---
    for idx = 1:length(activityCols)
        colIdx    = activityCols(idx);
        colHeader = varNames{colIdx};
        signalID  = colHeader;
        
        % Extract activity signal
        signal = dataTable{:, colIdx};
        if ~isnumeric(signal)
            signal = str2double(signal);
        end
        signal(~isfinite(signal)) = 0;
        
        % Append to group matrix
        signalsMatrix = [signalsMatrix, signal]; %#ok<AGROW>
        
        %% --- Wavelet Analysis for the Individual ---
        FB = cwtfilterbank('SignalLength', numel(signal), ...
                           'SamplingPeriod', minutes(samplingPeriod), ...
                           'PeriodLimits', [minutes(60), minutes(1590)], ...
                           'Wavelet', 'amor');
        [wt, periods, coi] = cwt(signal, 'FilterBank', FB);
        periods_hours = hours(periods);
        
        %% --- Plot Individual Scalogram ---
        fig = figure;
        pcolor(time_day, periods_hours, abs(wt));
        shading interp;
        colormap jet;
        colorbar;
        caxis auto;
        set(gca, 'box', 'off');
        set(gca, 'YTick', 0:4:26);
        hold on;
        % vertical lines at light-condition changes
        for k = 1:length(condChangeIdx)
            boundaryRow = condChangeIdx(k) + 1;
            x_line = time_day(boundaryRow);
            plot([x_line, x_line], [min(periods_hours) max(periods_hours)], ...
                 'w:', 'LineWidth', 1.5);
        end
        hold off;
        
        xlabel('Time (days)');
        ylabel('Period (hr)');
        title(sprintf('Scalogram - %s | %s', signalID, groupName));
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', ...
                 'FontName', 'Times New Roman');
        
        safeID     = strrep(signalID, ' ', '_');
        outFileInd = fullfile(indivFolder, ...
                        sprintf('Scalogram_%s_%s.jpg', groupName, safeID));
        print(fig, outFileInd, '-djpeg', '-r600');
        close(fig);
    end  % end individuals loop
    
    %% --- Group Average Scalogram (all individuals in this condition) ---
    if ~isempty(signalsMatrix)
        avg_signal = mean(signalsMatrix, 2);
        
        FB_avg = cwtfilterbank('SignalLength', numel(avg_signal), ...
                               'SamplingPeriod', minutes(samplingPeriod), ...
                               'PeriodLimits', [minutes(60), minutes(1590)], ...
                               'Wavelet', 'amor');
        [wt_avg, periods_avg, coi_avg] = cwt(avg_signal, 'FilterBank', FB_avg);
        periods_hours_avg = hours(periods_avg);
        
        fig_avg = figure;
        pcolor(time_day, periods_hours_avg, abs(wt_avg));
        shading interp;
        colormap jet;
        colorbar;
        caxis auto;
        set(gca, 'box', 'off');
        set(gca, 'YTick', 0:4:26);
        hold on;
        for k = 1:length(condChangeIdx)
            boundaryRow = condChangeIdx(k) + 1;
            x_line = time_day(boundaryRow);
            plot([x_line, x_line], [min(periods_hours_avg) max(periods_hours_avg)], ...
                 'w:', 'LineWidth', 1.5);
        end
        hold off;
        
        xlabel('Time (days)');
        ylabel('Period (hr)');
        title(sprintf('Scalogram - Average Signal | %s', groupName));
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', ...
                 'FontName', 'Times New Roman');
        
        outFileAvg = fullfile(groupFolder, ...
                        sprintf('Scalogram_Average_%s.jpg', groupName));
        print(fig_avg, outFileAvg, '-djpeg', '-r600');
        close(fig_avg);
    end
    
    fprintf('Completed analysis for condition group "%s".\n', groupName);
end

disp('Wavelet analysis completed. All individual and condition-specific average scalograms have been saved.');
