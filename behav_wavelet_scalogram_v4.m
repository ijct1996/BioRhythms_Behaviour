%% ------------------------------------------------------------------------
% Script: Wavelet Analysis of Mice Behaviour Data (Scalograms with Conditions
%         and Sex Separation)
%
% PURPOSE:
% Reads an Excel file containing mice behavioural data. After selecting:
%   - A time column (in hours; converted for plotting),
%   - A light condition column (for drawing vertical markers), and
%   - A number of condition groups (each with a user-specified name and the 
%     selection of activity columns corresponding to individuals in that condition),
%
% the script computes continuous wavelet transform scalograms for:
%   (a) each individual (activity column) in the condition, and
%   (b) the condition-specific average signal across those individuals.
%
% In this updated version, the script also separates outputs by sex.
% It detects the sex from each individual’s name (assuming it ends with “-M”
% for males or “-F” for females) and creates separate subfolders for saving
% scalograms. The average scalograms are computed separately for each sex.
%
% Each individual’s scalogram and the average scalogram are saved as 600 DPI 
% JPEG files in a dedicated subfolder for the condition and sex.
%
% NOTE: Adjust conversion factors, filter parameters, or aesthetic settings as needed.
%% ------------------------------------------------------------------------

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
varNames = dataTable.Properties.VariableNames;

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
% Ask how many condition groups (e.g., one for each experimental condition)
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

% For each condition group, prompt for a descriptive name and the corresponding
% activity columns (each column represents an individual in that condition).
condGroupNames = cell(numGroups,1);
condGroupColIdx = cell(numGroups,1);

for grp = 1:numGroups
    % Get the name for the condition group (used in file naming)
    promptMsg = {sprintf('Enter a name for condition group %d:', grp)};
    grpNameCell = inputdlg(promptMsg, 'Condition Group Name', [1 50]);
    if isempty(grpNameCell)
        error('No name entered for condition group %d. Aborting.', grp);
    end
    condGroupNames{grp} = grpNameCell{1};
    
    % Let the user select one or more activity data columns for this condition.
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
% Convert time (in hours) for plotting purposes.
time_hours = dataTable{:, timeIdx};
time_min   = time_hours * 60;             % Convert hours to minutes.
time_day   = time_min / (60 * 24);          % Convert minutes to days.
samplingPeriod = mean(diff(time_min));      % Sampling period in minutes.

% Extract the global Light Condition Column (for drawing vertical markers).
lightConditionVector = dataTable{:, lightCondIdx};

% Compute indices where the light condition changes (used for plotting markers).
condChangeIdx = find(diff(lightConditionVector) ~= 0);

%% --- Process Each Condition Group ---
% For each condition group, process all selected individuals and compute both
% individual and condition-specific average scalograms separated by sex.
for grp = 1:numGroups
    groupName = condGroupNames{grp};
    activityCols = condGroupColIdx{grp};
    
    % Create a subfolder for this condition group.
    groupFolder = fullfile(resultsFolder, ['Condition_' groupName]);
    if ~exist(groupFolder, 'dir')
       mkdir(groupFolder);
    end
    
    % Initialize matrices to store individual signals for averaging by sex.
    signalsMatrixM = [];  % For males
    signalsMatrixF = [];  % For females
    
    %% --- Loop over Each Activity Column (Individual) in the Group ---
    for idx = 1:length(activityCols)
        colIdx = activityCols(idx);
        colHeader = varNames{colIdx};
        signalID = colHeader;  % Use the column header as the signal identifier.
        
        % Determine sex from the signal ID (assumes ending with '-M' for male,
        % '-F' for female). If not found, mark as 'Unknown'.
        if endsWith(signalID, '-M')
            sexStr = 'Male';
        elseif endsWith(signalID, '-F')
            sexStr = 'Female';
        else
            sexStr = 'Unknown';
            warning('Sex could not be determined for %s. Assigned as Unknown.', signalID);
        end
        
        % Create a subfolder within the group folder for the given sex.
        sexFolder = fullfile(groupFolder, sexStr);
        if ~exist(sexFolder, 'dir')
            mkdir(sexFolder);
        end
        
        % Extract the activity signal for this individual.
        signal = dataTable{:, colIdx};
        if ~isnumeric(signal)
            signal = str2double(signal);
        end
        signal(~isfinite(signal)) = 0;
        
        % Append the signal to the appropriate sex-specific matrix.
        if strcmp(sexStr, 'Male')
            signalsMatrixM = [signalsMatrixM, signal];
        elseif strcmp(sexStr, 'Female')
            signalsMatrixF = [signalsMatrixF, signal];
        end
        
        %% --- Wavelet Analysis for the Individual's Signal ---
        % Set up a continuous wavelet filter bank with period limits
        % from 60 to 1590 minutes (approx. 1 hr to ~26.5 hrs) using the 'amor' wavelet.
        FB = cwtfilterbank('SignalLength', numel(signal), ...
                           'SamplingPeriod', minutes(samplingPeriod), ...
                           'PeriodLimits', [minutes(60), minutes(1590)], ...
                           'Wavelet', 'amor');
        [wt, periods, coi] = cwt(signal, 'FilterBank', FB);
        % Convert the period axis to hours for display.
        periods_hours = hours(periods);
        
        %% --- Plot the Individual's Magnitude Scalogram ---
        fig = figure;
        pcolor(time_day, periods_hours, abs(wt));
        shading interp;
        colormap jet;
        colorbar;
        caxis auto;  % Autoscale the colour limits.
        set(gca, 'box', 'off');
        set(gca, 'YTick', 0:4:26);
        hold on;
        % Draw vertical dotted white lines where the light condition changes.
        for k = 1:length(condChangeIdx)
            boundaryRow = condChangeIdx(k) + 1;
            x_line = time_day(boundaryRow);
            plot([x_line, x_line], [min(periods_hours) max(periods_hours)], 'w:', 'LineWidth', 1.5);
        end
        hold off;
        
        xlabel('Time (days)');
        ylabel('Period (hr)');
        title(sprintf('Scalogram - %s | %s | %s', signalID, groupName, sexStr));
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', ...
            'FontName', 'Times New Roman');
        
        % Save the individual scalogram into the appropriate sex folder.
        safeID = strrep(signalID, ' ', '_');
        outFileName = fullfile(sexFolder, sprintf('Scalogram_%s_%s_%s.jpg', groupName, safeID, sexStr));
        print(fig, outFileName, '-djpeg', '-r600');
        close(fig);
    end  % End individual loop
    
    %% --- Compute and Plot the Condition-Specific Average Scalograms by Sex ---
    % For Males
    if ~isempty(signalsMatrixM)
        avg_signal = mean(signalsMatrixM, 2);
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
            plot([x_line, x_line], [min(periods_hours_avg) max(periods_hours_avg)], 'w:', 'LineWidth', 1.5);
        end
        hold off;
        
        xlabel('Time (days)');
        ylabel('Period (hr)');
        title(sprintf('Scalogram - Average Signal | %s | Male', groupName));
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', ...
            'FontName', 'Times New Roman');
        
        maleFolder = fullfile(groupFolder, 'Male');
        if ~exist(maleFolder, 'dir')
            mkdir(maleFolder);
        end
        outFileName_avg = fullfile(maleFolder, sprintf('Scalogram_Average_%s_Male.jpg', groupName));
        print(fig_avg, outFileName_avg, '-djpeg', '-r600');
        close(fig_avg);
    end
    
    % For Females
    if ~isempty(signalsMatrixF)
        avg_signal = mean(signalsMatrixF, 2);
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
            plot([x_line, x_line], [min(periods_hours_avg) max(periods_hours_avg)], 'w:', 'LineWidth', 1.5);
        end
        hold off;
        
        xlabel('Time (days)');
        ylabel('Period (hr)');
        title(sprintf('Scalogram - Average Signal | %s | Female', groupName));
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', ...
            'FontName', 'Times New Roman');
        
        femaleFolder = fullfile(groupFolder, 'Female');
        if ~exist(femaleFolder, 'dir')
            mkdir(femaleFolder);
        end
        outFileName_avg = fullfile(femaleFolder, sprintf('Scalogram_Average_%s_Female.jpg', groupName));
        print(fig_avg, outFileName_avg, '-djpeg', '-r600');
        close(fig_avg);
    end
    
    fprintf('Completed analysis for condition group "%s".\n', groupName);
end

disp('Wavelet analysis completed. All individual and condition-specific average scalograms are saved by sex in the selected folder.');