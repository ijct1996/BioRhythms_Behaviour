%% ------------------------------------------------------------------------
% Script: Wavelet Analysis of Mice Behaviour Data with Power Spectra &
%         Peak Analysis Output to Excel (Condition Group Specific)
%
% PURPOSE:
% This script performs wavelet analysis on an Excel file containing mice
% behavioural data. The user is prompted to select:
%   - The time column (in hours, converted to minutes for analysis and to days 
%     for plotting on the x-axis)
%   - The light condition column (numeric values; used for drawing vertical markers)
%   - A number of condition groups. For each condition group, the user specifies a
%     name and selects one or more activity data columns (each representing an individual
%     within that condition).
%
% For each individual signal within a condition group, the script:
%   - Computes the continuous wavelet transform using a filter bank covering periods
%     from 60 minutes (1 hr) to 1590 minutes (~26.5 hr) with the 'amor' wavelet.
%   - Computes a global average power spectrum (mean of log10(power) over time) and plots
%     it.
%   - For each unique light condition, computes a condition-specific average power 
%     spectrum, detects peaks using findpeaks (with a minimum peak prominence of 0),
%     sorts the peaks (strongest first), plots the spectrum, and saves the plot.
%
% For each condition group, the script also computes the group average signal across all
% individuals and outputs the corresponding power spectra (global and condition specific,
% with peak detection).
%
% Peak summary data (SignalID, PeakRank, PeakPeriod_hr, PeakValue_log10, PeakProminence, 
% PeakWidth) for each light condition (across all signals, including individual and group 
% averages) is accumulated and written to an Excel file (one sheet per condition).
%
% Aesthetic settings include: no grid, tick direction out, Times New Roman font, and no
% top/right box.
%
% CHANGELOG:
% - 2025-11-09: Original version included sex detection based on '-M'/'-F' suffixes and
%               sex specific folders and averages (Male, Female, Unknown).
% - 2025-11-09: Removed sex detection and sex specific folder structure. All selected
%               individuals in a condition group are treated as one pool; a single
%               group average is computed. Simplified folder structure.
% - 2025-11-09: Fixed globalPowerMax calculation so that the margin is added once after
%               scanning all spectra instead of incrementing inside the loop.
%% ------------------------------------------------------------------------

%% --- Clear Environment ---
clearvars; close all; clc;

%% --- Select the Input Excel File ---
[fileName, filePath] = uigetfile('*.xlsx', 'Select the source Excel file');
if isequal(fileName, 0)
    error('No file selected. Script aborted.');
end
sourceFile = fullfile(filePath, fileName);

%% --- Read the Source Data (Preserving Original Column Headers) ---
dataTable = readtable(sourceFile, 'VariableNamingRule', 'preserve');
varNames  = dataTable.Properties.VariableNames;

%% --- User Selections: Time & Light Condition Columns ---
% Select the Time Column (assumed to be in hours)
[timeIdx, tf] = listdlg('PromptString', 'Select the Time Column (in hours)', ...
                        'SelectionMode', 'single', 'ListString', varNames);
if ~tf
    error('Time column not selected. Aborting.');
end

% Select the Light Condition Column (numeric values; used for drawing vertical markers)
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
numGroups = str2double(strtrim(numCondStr{1}));
if isnan(numGroups) || numGroups < 1 || mod(numGroups,1) ~= 0
    error('Invalid number of condition groups entered. Aborting.');
end

% For each condition group, prompt for a name and let the user select one or more
% activity data columns.
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

%% --- Extract and Convert Time Vector ---
% Time is assumed to be in hours.
time_hours      = dataTable{:, timeIdx};
time_min        = time_hours * 60;               % Convert hours to minutes.
time_day        = time_min / (60 * 24);          % For plotting: convert minutes to days.
samplingPeriod  = mean(diff(time_min));          % Sampling period (in minutes).

%% --- Extract Light Condition Column ---
conditionVector = dataTable{:, lightCondIdx};

%% --- Set Global y-Axis Limit Dynamically ---
globalYMax = hours(minutes(1590));   % ~26.5 hr

%% --- Determine Global Maximum Power (for x-axis) ---
allActivityIndices = unique([condGroupColIdx{:}]);
globalPowerMax     = -inf;
uniqueCond         = unique(conditionVector);
nCond              = numel(uniqueCond);

for idx = allActivityIndices'
    signal = dataTable{:, idx};
    if ~isvector(signal)
        signal = signal(:);
    end
    if ~isnumeric(signal)
        signal = str2double(signal);
    end
    signal(~isfinite(signal)) = 0;
    
    FB = cwtfilterbank('SignalLength', numel(signal), ...
                       'SamplingPeriod', minutes(samplingPeriod), ...
                       'PeriodLimits', [minutes(60), minutes(1590)], ...
                       'Wavelet', 'amor');
    [wt, periods, coi] = cwt(signal, 'FilterBank', FB);
    powerSpec = abs(wt).^2;
    
    % Global average power spectrum over time.
    avgPowerSpectrum = mean(log10(powerSpec), 2);
    globalPowerMax   = max(globalPowerMax, max(avgPowerSpectrum));
    
    % Condition specific averages.
    for j = 1:nCond
        condIdx = conditionVector == uniqueCond(j);
        if any(condIdx)
            condAvgPowerSpectrum = mean(log10(powerSpec(:, condIdx)), 2);
            globalPowerMax       = max(globalPowerMax, max(condAvgPowerSpectrum));
        end
    end
end

% Add margin once and ensure lower bound is at least 0.
if ~isfinite(globalPowerMax)
    globalPowerMax = 0;
else
    globalPowerMax = globalPowerMax + 0.8;
    if globalPowerMax < 0
        globalPowerMax = 0;
    end
end

%% --- Prepare for Peak Accumulation ---
% Each cell holds: {SignalID, PeakRank, PeakPeriod_hr, PeakValue_log10, PeakProminence, PeakWidth}
peakResults = cell(nCond, 1);
for i = 1:nCond
    peakResults{i} = {};
end

%% --- Process Each Condition Group ---
for grp = 1:numGroups
    groupName    = condGroupNames{grp};
    activityCols = condGroupColIdx{grp};
    
    % Folder for this group
    groupFolder = fullfile(resultsFolder, ['ConditionGroup_' groupName]);
    if ~exist(groupFolder, 'dir')
        mkdir(groupFolder);
    end
    
    % Optional subfolder for individual spectra
    indivFolder = fullfile(groupFolder, 'Individuals');
    if ~exist(indivFolder, 'dir')
        mkdir(indivFolder);
    end
    
    % Matrix for group average (all individuals)
    groupSignalsMatrix = [];
    
    % Process each individual
    for i = 1:length(activityCols)
        colIdx    = activityCols(i);
        colHeader = varNames{colIdx};
        signalID  = sprintf('%s_%s', groupName, colHeader);
        
        % Extract the signal
        signal = dataTable{:, colIdx};
        if ~isnumeric(signal)
            signal = str2double(signal);
        end
        signal(~isfinite(signal)) = 0;
        
        % Store for group average
        groupSignalsMatrix = [groupSignalsMatrix, signal]; %#ok<AGROW>
        
        %% --- Wavelet Analysis for Individual Signal ---
        FB = cwtfilterbank('SignalLength', numel(signal), ...
                           'SamplingPeriod', minutes(samplingPeriod), ...
                           'PeriodLimits', [minutes(60), minutes(1590)], ...
                           'Wavelet', 'amor');
        [wt, periods, coi] = cwt(signal, 'FilterBank', FB);
        powerSpec          = abs(wt).^2;
        avgPowerSpectrum   = mean(log10(powerSpec), 2);
        periods_hours      = hours(periods);
        
        %% --- Global Average Power Spectrum Plot (Individual) ---
        figGlobal = figure;
        plot(avgPowerSpectrum, periods_hours, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})');
        ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - %s', signalID), 'Interpreter', 'none');
        set(gca, 'XGrid','off','YGrid','off','TickDir','out', ...
                 'FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]);
        ylim([0, globalYMax]);
        
        safeID = strrep(signalID, ' ', '_');
        outFileNameGlobal = fullfile(indivFolder, ...
            sprintf('PowerSpectrum_%s_Global.jpg', safeID));
        print(figGlobal, outFileNameGlobal, '-djpeg', '-r600');
        close(figGlobal);
        
        %% --- Condition Specific Spectra & Peaks (Individual) ---
        for j = 1:nCond
            currCond   = uniqueCond(j);
            condIdx    = conditionVector == currCond;
            if ~any(condIdx)
                continue;
            end
            
            condAvgPowerSpectrum = mean(log10(powerSpec(:, condIdx)), 2);
            [condPeaks, condLocs, condWidths, condProms] = ...
                findpeaks(condAvgPowerSpectrum, 'MinPeakProminence', 0);
            if isempty(condPeaks)
                continue;
            end
            
            [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
            sortedLocs             = condLocs(sortIdx);
            sortedWidths           = condWidths(sortIdx);
            sortedProms            = condProms(sortIdx);
            sortedPeriods          = periods_hours(sortedLocs);
            
            figCond = figure;
            plot(condAvgPowerSpectrum, periods_hours, '-k', 'LineWidth', 1.5);
            xlabel('Power (log_{10})');
            ylabel('Period (hr)');
            title(sprintf('Power Spectrum - %s, Light Cond: %d', signalID, currCond), ...
                'Interpreter', 'none');
            set(gca, 'XGrid','off','YGrid','off','TickDir','out', ...
                     'FontName','Times New Roman','box','off');
            xlim([0, globalPowerMax]);
            ylim([0, globalYMax]);
            
            outFileNameCond = fullfile(indivFolder, ...
                sprintf('PowerSpectrum_%s_Condition_%d.jpg', safeID, currCond));
            print(figCond, outFileNameCond, '-djpeg', '-r600');
            close(figCond);
            
            nPeaks = numel(sortedPeaks);
            for p = 1:nPeaks
                newRow = {signalID, p, sortedPeriods(p), sortedPeaks(p), ...
                          sortedProms(p), sortedWidths(p)};
                peakResults{j} = [peakResults{j}; newRow];
            end
        end
    end % individuals loop
    
    %% --- Group Average Across All Individuals (this condition group) ---
    if ~isempty(groupSignalsMatrix)
        groupAvgSignal = mean(groupSignalsMatrix, 2);
        
        FB_avg = cwtfilterbank('SignalLength', numel(groupAvgSignal), ...
                               'SamplingPeriod', minutes(samplingPeriod), ...
                               'PeriodLimits', [minutes(60), minutes(1590)], ...
                               'Wavelet', 'amor');
        [wt_avg, periods_avg, coi_avg] = cwt(groupAvgSignal, 'FilterBank', FB_avg);
        powerSpec_avg        = abs(wt_avg).^2;
        avgPowerSpectrum_avg = mean(log10(powerSpec_avg), 2);
        periods_hours_avg    = hours(periods_avg);
        
        % Global average spectrum for group
        figAvgGlobal = figure;
        plot(avgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})');
        ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - Average for %s', groupName), ...
            'Interpreter', 'none');
        set(gca, 'XGrid','off','YGrid','off','TickDir','out', ...
                 'FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]);
        ylim([0, globalYMax]);
        
        safeGroup = strrep(groupName, ' ', '_');
        outFileNameGroupGlobal = fullfile(groupFolder, ...
            sprintf('PowerSpectrum_%s_Global_Average.jpg', safeGroup));
        print(figAvgGlobal, outFileNameGroupGlobal, '-djpeg', '-r600');
        close(figAvgGlobal);
        
        % Condition specific spectra + peaks for group average
        for j = 1:nCond
            currCond = uniqueCond(j);
            condIdx  = conditionVector == currCond;
            if ~any(condIdx)
                continue;
            end
            
            condAvgPowerSpectrum_avg = mean(log10(powerSpec_avg(:, condIdx)), 2);
            [condPeaks, condLocs, condWidths, condProms] = ...
                findpeaks(condAvgPowerSpectrum_avg, 'MinPeakProminence', 0);
            if isempty(condPeaks)
                continue;
            end
            
            [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
            sortedLocs             = condLocs(sortIdx);
            sortedWidths           = condWidths(sortIdx);
            sortedProms            = condProms(sortIdx);
            sortedPeriods          = periods_hours_avg(sortedLocs);
            
            figCondAvg = figure;
            plot(condAvgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
            xlabel('Power (log_{10})');
            ylabel('Period (hr)');
            title(sprintf('Power Spectrum - Avg for %s, Light Cond: %d', ...
                groupName, currCond), 'Interpreter', 'none');
            set(gca, 'XGrid','off','YGrid','off','TickDir','out', ...
                     'FontName','Times New Roman','box','off');
            xlim([0, globalPowerMax]);
            ylim([0, globalYMax]);
            
            outFileNameCondAvg = fullfile(groupFolder, ...
                sprintf('PowerSpectrum_%s_Average_Condition_%d.jpg', ...
                    safeGroup, currCond));
            print(figCondAvg, outFileNameCondAvg, '-djpeg', '-r600');
            close(figCondAvg);
            
            nPeaks = numel(sortedPeaks);
            avgID  = [groupName '_Average'];
            for p = 1:nPeaks
                newRow = {avgID, p, sortedPeriods(p), sortedPeaks(p), ...
                          sortedProms(p), sortedWidths(p)};
                peakResults{j} = [peakResults{j}; newRow];
            end
        end
    end
    
    fprintf('Completed analysis for condition group "%s".\n', groupName);
end

%% --- Write Excel File for Peak Summaries ---
excelFileName = fullfile(resultsFolder, 'PeakResults.xlsx');
headers       = {'SignalID', 'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', ...
                 'PeakProminence', 'PeakWidth'};

for j = 1:nCond
    if isempty(peakResults{j})
        continue;
    end
    T = cell2table(peakResults{j}, 'VariableNames', headers);
    T = sortrows(T, 'PeakValue_log10', 'descend');
    sheetName = sprintf('Condition_%d', uniqueCond(j));
    writetable(T, excelFileName, 'Sheet', sheetName);
end

disp('Wavelet analysis completed. Global and condition specific power spectra (with peak analyses) have been saved as JPEG images, and peak results have been written to an Excel file.');
