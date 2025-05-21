%% ------------------------------------------------------------------------
% Script: Wavelet Analysis of Mice Behaviour Data with Power Spectra &
%         Peak Analysis Output to Excel (Condition Group Specific) with Sex
%         Separation
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
% The analysis is performed separately for males and females. The sex is determined
% by examining the column header (expected to end with '-M' for males or '-F' for females).
%
% For each condition group, the script also computes the group average signal across all
% individuals separated by sex and outputs the corresponding power spectra (global and 
% condition-specific, with peak detection).
%
% Peak summary data (SignalID, PeakRank, PeakPeriod_hr, PeakValue_log10, PeakProminence, 
% PeakWidth) for each light condition (across all signals, including individual and group 
% averages) is accumulated and written to an Excel file (one sheet per condition).
%
% Aesthetic settings include: no grid, tick direction out, Times New Roman font, and no
% top/right box.
%
% NOTE: Adjust conversion factors, filter parameters, or aesthetic settings as needed.
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
varNames = dataTable.Properties.VariableNames;

%% --- User Selections: Time & Light Condition Columns ---
% Select the Time Column (assumed to be in hours)
[timeIdx, tf] = listdlg('PromptString', 'Select the Time Column (in hours)', ...
                        'SelectionMode', 'single', 'ListString', varNames);
if ~tf, error('Time column not selected. Aborting.'); end

% Select the Light Condition Column (numeric values; used for drawing vertical markers)
[lightCondIdx, tf] = listdlg('PromptString', 'Select the Light Condition Column', ...
                             'SelectionMode', 'single', 'ListString', varNames);
if ~tf, error('Light condition column not selected. Aborting.'); end

%% --- Prompt for Condition Groups ---
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

% For each condition group, prompt for a name and let the user select one or more
% activity data columns.
condGroupNames = cell(numGroups,1);
condGroupColIdx = cell(numGroups,1);
for grp = 1:numGroups
    promptMsg = {sprintf('Enter a name for condition group %d:', grp)};
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
if resultsFolder == 0, error('No output folder selected. Script aborted.'); end

%% --- Extract and Convert Time Vector ---
% Time is assumed to be in hours.
time_hours = dataTable{:, timeIdx};
time_min   = time_hours * 60;                % Convert hours to minutes.
time_day   = time_min / (60*24);               % For plotting: convert minutes to days.
samplingPeriod = mean(diff(time_min));         % Sampling period (in minutes).

%% --- Extract Light Condition Column ---
conditionVector = dataTable{:, lightCondIdx};

%% --- Set Global y-Axis Limit Dynamically ---
% Use the upper period limit of the wavelet analysis filter bank.
globalYMax = hours(minutes(1590));   % This equals approximately 26.5 hr.

%% --- Determine Global Maximum Power ---
% This loop processes every selected signal, ensures it is a vector, and
% computes the maximum value among all global and condition-specific average
% power spectra (on the log10 scale).
allActivityIndices = unique([condGroupColIdx{:}]);
globalPowerMax = -inf;
uniqueCond = unique(conditionVector);
nCond = numel(uniqueCond);
for idx = allActivityIndices'
    signal = dataTable{:, idx};
    if ~isvector(signal)
        signal = signal(:);  % force as column vector
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
    globalPowerMax_0 = max(globalPowerMax, max(avgPowerSpectrum));
    globalPowerMax = globalPowerMax_0 + 0.8; 
    % Also check condition-specific averages.
    for j = 1:nCond
         condIndices = find(conditionVector == uniqueCond(j));
         if ~isempty(condIndices)
             condAvgPowerSpectrum = mean(log10(powerSpec(:, condIndices)), 2);
             globalPowerMax = max(globalPowerMax, max(condAvgPowerSpectrum));
         end
    end
end
% Ensure that the lower bound is at least 0.
if globalPowerMax < 0
    globalPowerMax = 0;
end

%% --- Prepare for Peak Accumulation ---
% Each cell will hold peak summary data: {SignalID, PeakRank, PeakPeriod_hr, PeakValue_log10, PeakProminence, PeakWidth}.
peakResults = cell(nCond, 1);
for i = 1:nCond
    peakResults{i} = {};
end

%% --- Process Each Condition Group ---
for grp = 1:numGroups
    groupName = condGroupNames{grp};
    activityCols = condGroupColIdx{grp};  % Activity columns for this group.
    
    % Create a dedicated subfolder for this group.
    groupFolder = fullfile(resultsFolder, ['ConditionGroup_' groupName]);
    if ~exist(groupFolder, 'dir')
       mkdir(groupFolder);
    end
    % Create subfolders for each sex.
    maleFolder = fullfile(groupFolder, 'Male');
    if ~exist(maleFolder, 'dir'), mkdir(maleFolder); end
    femaleFolder = fullfile(groupFolder, 'Female');
    if ~exist(femaleFolder, 'dir'), mkdir(femaleFolder); end
    unknownFolder = fullfile(groupFolder, 'Unknown');
    if ~exist(unknownFolder, 'dir'), mkdir(unknownFolder); end

    % Initialize matrices for group averages by sex.
    groupSignalsMatrixM = [];
    groupSignalsMatrixF = [];
    groupSignalsMatrixUnknown = [];
    
    % Process each individual.
    for i = 1:length(activityCols)
        colIdx = activityCols(i);
        colHeader = varNames{colIdx};
        
        % Determine sex.
        if endsWith(colHeader, '-M')
            sexStr = 'Male';
        elseif endsWith(colHeader, '-F')
            sexStr = 'Female';
        else
            sexStr = 'Unknown';
            warning('Sex could not be determined for %s. Assigned as Unknown.', colHeader);
        end
        
        % Construct signal ID.
        signalID = sprintf('%s_%s_%s', groupName, colHeader, sexStr);
        
        % Determine output folder based on sex.
        switch sexStr
            case 'Male'
                outFolder = maleFolder;
            case 'Female'
                outFolder = femaleFolder;
            otherwise
                outFolder = unknownFolder;
        end
        
        % Extract the signal.
        signal = dataTable{:, colIdx};
        if ~isnumeric(signal)
            signal = str2double(signal);
        end
        signal(~isfinite(signal)) = 0;
        
        % Append signal for group average.
        switch sexStr
            case 'Male'
                groupSignalsMatrixM = [groupSignalsMatrixM, signal];
            case 'Female'
                groupSignalsMatrixF = [groupSignalsMatrixF, signal];
            otherwise
                groupSignalsMatrixUnknown = [groupSignalsMatrixUnknown, signal];
        end
        
        %% --- Wavelet Analysis for Individual Signal ---
        FB = cwtfilterbank('SignalLength', numel(signal), ...
                           'SamplingPeriod', minutes(samplingPeriod), ...
                           'PeriodLimits', [minutes(60), minutes(1590)], ...
                           'Wavelet', 'amor');
        [wt, periods, coi] = cwt(signal, 'FilterBank', FB);
        powerSpec = abs(wt).^2;
        avgPowerSpectrum = mean(log10(powerSpec), 2);
        periods_hours = hours(periods);
        
        %% --- Global Average Power Spectrum Plot for Individual ---
        figGlobal = figure;
        plot(avgPowerSpectrum, periods_hours, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})');
        ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - %s', signalID), 'Interpreter', 'none');
        set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]);
        ylim([0, globalYMax]);
        safeID = strrep(signalID, ' ', '_');
        outFileNameGlobal = fullfile(outFolder, sprintf('PowerSpectrum_%s_Global.jpg', safeID));
        print(figGlobal, outFileNameGlobal, '-djpeg', '-r600');
        close(figGlobal);
        
        %% --- Condition-Specific Power Spectrum & Peak Detection for Individual ---
        for j = 1:nCond
            currCond = uniqueCond(j);
            condIndices = find(conditionVector == currCond);
            if ~isempty(condIndices)
                condAvgPowerSpectrum = mean(log10(powerSpec(:, condIndices)), 2);
                [condPeaks, condPeakLocs, condWidths, condProminences] = findpeaks(condAvgPowerSpectrum, 'MinPeakProminence', 0);
                if isempty(condPeaks)
                    continue;
                end
                [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
                sortedPeakLocs = condPeakLocs(sortIdx);
                sortedWidths = condWidths(sortIdx);
                sortedProminences = condProminences(sortIdx);
                sortedPeakPeriods = periods_hours(sortedPeakLocs);
                
                figCond = figure;
                plot(condAvgPowerSpectrum, periods_hours, '-k', 'LineWidth', 1.5);
                xlabel('Power (log_{10})');
                ylabel('Period (hr)');
                title(sprintf('Power Spectrum - %s, Light Cond: %d', signalID, currCond), 'Interpreter', 'none');
                set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
                xlim([0, globalPowerMax]);
                ylim([0, globalYMax]);
                outFileNameCond = fullfile(outFolder, sprintf('PowerSpectrum_%s_Condition_%d.jpg', safeID, currCond));
                print(figCond, outFileNameCond, '-djpeg', '-r600');
                close(figCond);
                
                nPeaks = numel(sortedPeaks);
                for p = 1:nPeaks
                    newRow = {signalID, p, sortedPeakPeriods(p), sortedPeaks(p), sortedProminences(p), sortedWidths(p)};
                    peakResults{j} = [peakResults{j}; newRow];
                end
            end
        end
    end % End individuals loop.
    
    %% --- Process Group Average (Across Individuals) by Sex ---
    % For Males:
    if ~isempty(groupSignalsMatrixM)
        groupAvgSignal = mean(groupSignalsMatrixM, 2);
        FB_avg = cwtfilterbank('SignalLength', numel(groupAvgSignal), ...
                               'SamplingPeriod', minutes(samplingPeriod), ...
                               'PeriodLimits', [minutes(60), minutes(1590)], ...
                               'Wavelet', 'amor');
        [wt_avg, periods_avg, coi_avg] = cwt(groupAvgSignal, 'FilterBank', FB_avg);
        powerSpec_avg = abs(wt_avg).^2;
        avgPowerSpectrum_avg = mean(log10(powerSpec_avg), 2);
        periods_hours_avg = hours(periods_avg);
        
        figAvgGlobal = figure;
        plot(avgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})');
        ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - Average for %s_Male', groupName), 'Interpreter', 'none');
        set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]);
        ylim([0, globalYMax]);
        safeGroup = strrep(groupName, ' ', '_');
        outFileNameGroupGlobal = fullfile(maleFolder, sprintf('PowerSpectrum_%s_Global_Average_Male.jpg', safeGroup));
        print(figAvgGlobal, outFileNameGroupGlobal, '-djpeg', '-r600');
        close(figAvgGlobal);
        
        for j = 1:nCond
            currCond = uniqueCond(j);
            condIndices = find(conditionVector == currCond);
            if ~isempty(condIndices)
                condAvgPowerSpectrum_avg = mean(log10(powerSpec_avg(:, condIndices)), 2);
                [condPeaks, condPeakLocs, condWidths, condProminences] = findpeaks(condAvgPowerSpectrum_avg, 'MinPeakProminence', 0);
                if isempty(condPeaks)
                    continue;
                end
                [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
                sortedPeakLocs = condPeakLocs(sortIdx);
                sortedWidths = condWidths(sortIdx);
                sortedProminences = condProminences(sortIdx);
                sortedPeakPeriods = periods_hours_avg(sortedPeakLocs);
                
                figCondAvg = figure;
                plot(condAvgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
                xlabel('Power (log_{10})');
                ylabel('Period (hr)');
                title(sprintf('Power Spectrum - Avg for %s_Male, Light Cond: %d', groupName, currCond), 'Interpreter', 'none');
                set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
                xlim([0, globalPowerMax]);
                ylim([0, globalYMax]);
                outFileNameCondAvg = fullfile(maleFolder, sprintf('PowerSpectrum_%s_Average_Condition_%d_Male.jpg', safeGroup, currCond));
                print(figCondAvg, outFileNameCondAvg, '-djpeg', '-r600');
                close(figCondAvg);
                
                nPeaks = numel(sortedPeaks);
                for p = 1:nPeaks
                    newRow = { [groupName '_Average_Male'], p, sortedPeakPeriods(p), sortedPeaks(p), sortedProminences(p), sortedWidths(p) };
                    peakResults{j} = [peakResults{j}; newRow];
                end
            end
        end
    end
    
    % For Females:
    if ~isempty(groupSignalsMatrixF)
        groupAvgSignal = mean(groupSignalsMatrixF, 2);
        FB_avg = cwtfilterbank('SignalLength', numel(groupAvgSignal), ...
                               'SamplingPeriod', minutes(samplingPeriod), ...
                               'PeriodLimits', [minutes(60), minutes(1590)], ...
                               'Wavelet', 'amor');
        [wt_avg, periods_avg, coi_avg] = cwt(groupAvgSignal, 'FilterBank', FB_avg);
        powerSpec_avg = abs(wt_avg).^2;
        avgPowerSpectrum_avg = mean(log10(powerSpec_avg), 2);
        periods_hours_avg = hours(periods_avg);
        
        figAvgGlobal = figure;
        plot(avgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})');
        ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - Average for %s_Female', groupName), 'Interpreter', 'none');
        set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]);
        ylim([0, globalYMax]);
        safeGroup = strrep(groupName, ' ', '_');
        outFileNameGroupGlobal = fullfile(femaleFolder, sprintf('PowerSpectrum_%s_Global_Average_Female.jpg', safeGroup));
        print(figAvgGlobal, outFileNameGroupGlobal, '-djpeg', '-r600');
        close(figAvgGlobal);
        
        for j = 1:nCond
            currCond = uniqueCond(j);
            condIndices = find(conditionVector == currCond);
            if ~isempty(condIndices)
                condAvgPowerSpectrum_avg = mean(log10(powerSpec_avg(:, condIndices)), 2);
                [condPeaks, condPeakLocs, condWidths, condProminences] = findpeaks(condAvgPowerSpectrum_avg, 'MinPeakProminence', 0);
                if isempty(condPeaks)
                    continue;
                end
                [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
                sortedPeakLocs = condPeakLocs(sortIdx);
                sortedWidths = condWidths(sortIdx);
                sortedProminences = condProminences(sortIdx);
                sortedPeakPeriods = periods_hours_avg(sortedPeakLocs);
                
                figCondAvg = figure;
                plot(condAvgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
                xlabel('Power (log_{10})');
                ylabel('Period (hr)');
                title(sprintf('Power Spectrum - Avg for %s_Female, Light Cond: %d', groupName, currCond), 'Interpreter', 'none');
                set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
                xlim([0, globalPowerMax]);
                ylim([0, globalYMax]);
                outFileNameCondAvg = fullfile(femaleFolder, sprintf('PowerSpectrum_%s_Average_Condition_%d_Female.jpg', safeGroup, currCond));
                print(figCondAvg, outFileNameCondAvg, '-djpeg', '-r600');
                close(figCondAvg);
                
                nPeaks = numel(sortedPeaks);
                for p = 1:nPeaks
                    newRow = { [groupName '_Average_Female'], p, sortedPeakPeriods(p), sortedPeaks(p), sortedProminences(p), sortedWidths(p) };
                    peakResults{j} = [peakResults{j}; newRow];
                end
            end
        end
    end
    
    % For Unknown sex:
    if ~isempty(groupSignalsMatrixUnknown)
        groupAvgSignal = mean(groupSignalsMatrixUnknown, 2);
        FB_avg = cwtfilterbank('SignalLength', numel(groupAvgSignal), ...
                               'SamplingPeriod', minutes(samplingPeriod), ...
                               'PeriodLimits', [minutes(60), minutes(1590)], ...
                               'Wavelet', 'amor');
        [wt_avg, periods_avg, coi_avg] = cwt(groupAvgSignal, 'FilterBank', FB_avg);
        powerSpec_avg = abs(wt_avg).^2;
        avgPowerSpectrum_avg = mean(log10(powerSpec_avg), 2);
        periods_hours_avg = hours(periods_avg);
        
        figAvgGlobal = figure;
        plot(avgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
        xlabel('Power (log_{10})');
        ylabel('Period (hr)');
        title(sprintf('Global Power Spectrum - Average for %s_Unknown', groupName), 'Interpreter', 'none');
        set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
        xlim([0, globalPowerMax]);
        ylim([0, globalYMax]);
        safeGroup = strrep(groupName, ' ', '_');
        outFileNameGroupGlobal = fullfile(unknownFolder, sprintf('PowerSpectrum_%s_Global_Average_Unknown.jpg', safeGroup));
        print(figAvgGlobal, outFileNameGroupGlobal, '-djpeg', '-r600');
        close(figAvgGlobal);
        
        for j = 1:nCond
            currCond = uniqueCond(j);
            condIndices = find(conditionVector == currCond);
            if ~isempty(condIndices)
                condAvgPowerSpectrum_avg = mean(log10(powerSpec_avg(:, condIndices)), 2);
                [condPeaks, condPeakLocs, condWidths, condProminences] = findpeaks(condAvgPowerSpectrum_avg, 'MinPeakProminence', 0);
                if isempty(condPeaks)
                    continue;
                end
                [sortedPeaks, sortIdx] = sort(condPeaks, 'descend');
                sortedPeakLocs = condPeakLocs(sortIdx);
                sortedWidths = condWidths(sortIdx);
                sortedProminences = condProminences(sortIdx);
                sortedPeakPeriods = periods_hours_avg(sortedPeakLocs);
                
                figCondAvg = figure;
                plot(condAvgPowerSpectrum_avg, periods_hours_avg, '-k', 'LineWidth', 1.5);
                xlabel('Power (log_{10})');
                ylabel('Period (hr)');
                title(sprintf('Power Spectrum - Avg for %s_Unknown, Light Cond: %d', groupName, currCond), 'Interpreter', 'none');
                set(gca, 'XGrid','off','YGrid','off','TickDir','out','FontName','Times New Roman','box','off');
                xlim([0, globalPowerMax]);
                ylim([0, globalYMax]);
                outFileNameCondAvg = fullfile(unknownFolder, sprintf('PowerSpectrum_%s_Average_Condition_%d_Unknown.jpg', safeGroup, currCond));
                print(figCondAvg, outFileNameCondAvg, '-djpeg', '-r600');
                close(figCondAvg);
                
                nPeaks = numel(sortedPeaks);
                for p = 1:nPeaks
                    newRow = { [groupName '_Average_Unknown'], p, sortedPeakPeriods(p), sortedPeaks(p), sortedProminences(p), sortedWidths(p) };
                    peakResults{j} = [peakResults{j}; newRow];
                end
            end
        end
    end
    fprintf('Completed analysis for condition group "%s".\n', groupName);
end

%% --- Write Excel File for Peak Summaries ---
excelFileName = fullfile(resultsFolder, 'PeakResults.xlsx');
headers = {'SignalID', 'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', 'PeakProminence', 'PeakWidth'};
for j = 1:nCond
    if isempty(peakResults{j})
        continue;
    end
    T = cell2table(peakResults{j}, 'VariableNames', headers);
    T = sortrows(T, 'PeakValue_log10', 'descend');
    sheetName = sprintf('Condition_%d', uniqueCond(j));
    writetable(T, excelFileName, 'Sheet', sheetName);
end

disp('Wavelet analysis completed. Global and condition-specific power spectra (with peak analyses) have been saved as JPEG images, and peak results have been written to an Excel file.');
