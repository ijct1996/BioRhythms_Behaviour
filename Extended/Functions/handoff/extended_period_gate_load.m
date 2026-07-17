function [allCandidates, LoadSummary] = extended_period_gate_load(sumFiles, LOG)
%EXTENDED_PERIOD_GATE_LOAD Load PeriodCandidates_Long from WP_Summary__*.mat.

    if nargin < 2, LOG = -1; end

    allCandidates = table();
    loadRows = {};
    loadHdr = {'SummaryMat', 'Loaded', 'NRows', 'Message'};

    for i = 1:numel(sumFiles)
        matPath = fullfile(sumFiles(i).folder, sumFiles(i).name);
        try
            S = load(matPath);
            if isfield(S, 'pkgS') && isfield(S.pkgS, 'tables') && ...
                    isfield(S.pkgS.tables, 'PeriodCandidates_Long')
                T = S.pkgS.tables.PeriodCandidates_Long;
            elseif isfield(S, 'PeriodCandidates_Long')
                T = S.PeriodCandidates_Long;
            else
                error('No pkgS.tables.PeriodCandidates_Long table found.');
            end

            T = standardise_candidate_table_(T, matPath);
            allCandidates = vertcat_compatible_(allCandidates, T); %#ok<AGROW>
            loadRows(end+1, :) = {matPath, true, height(T), 'OK'}; %#ok<SAGROW>
            extended_period_gate_log(LOG, 'Loaded %s (%d rows).', matPath, height(T));
        catch ME
            loadRows(end+1, :) = {matPath, false, 0, ME.message}; %#ok<SAGROW>
            extended_period_gate_log(LOG, 'FAILED loading %s: %s', matPath, ME.message);
        end
    end

    LoadSummary = cell2table(loadRows, 'VariableNames', loadHdr);
end

%% ------------------------------------------------------------------------
function T = standardise_candidate_table_(T, sourcePackage)
    if ~istable(T)
        error('PeriodCandidates_Long object is not a table.');
    end

    strVars = {'File', 'SignalID', 'ConditionParsed', 'Source', 'HSubResidualMode', ...
        'Phase', 'BandName', 'CandidateID', 'QCReason'};
    numVars = {'Photoperiod_h', 'CandidateRank', 'MedianRidgePeriod_h', 'IQR_RidgePeriod_h', ...
        'MeanBandPower_log10', 'SDBandPower_log10', 'MeanRidgePower_log10', 'SDRidgePower_log10', ...
        'RidgeCoverageFrac', 'COIValidFrac', 'ValidPointCount', 'TotalPointCount'};
    logVars = {'PassQC'};

    for i = 1:numel(strVars)
        v = strVars{i};
        if ~ismember(v, T.Properties.VariableNames)
            T.(v) = strings(height(T), 1);
        else
            T.(v) = string(T.(v));
        end
    end

    for i = 1:numel(numVars)
        v = numVars{i};
        if ~ismember(v, T.Properties.VariableNames)
            T.(v) = nan(height(T), 1);
        else
            T.(v) = to_double_col_(T.(v));
        end
    end

    for i = 1:numel(logVars)
        v = logVars{i};
        if ~ismember(v, T.Properties.VariableNames)
            T.(v) = false(height(T), 1);
        else
            T.(v) = to_logical_col_(T.(v));
        end
    end

    missingCID = ismissing(T.CandidateID) | strlength(T.CandidateID) == 0;
    if any(missingCID)
        for r = find(missingCID(:))'
            T.CandidateID(r) = sprintf('%s|%s|%s|%s|%s|%s|%03d', ...
                char(T.File(r)), char(T.SignalID(r)), char(T.Source(r)), ...
                char(T.HSubResidualMode(r)), char(T.Phase(r)), char(T.BandName(r)), r);
        end
    end

    if ~ismember('SourcePackage', T.Properties.VariableNames)
        T.SourcePackage = repmat(string(sourcePackage), height(T), 1);
    else
        T.SourcePackage = string(T.SourcePackage);
    end

    keepVars = [{'SourcePackage'}, strVars, numVars, logVars];
    T = T(:, keepVars);
end

function out = vertcat_compatible_(A, B)
    if isempty(A)
        out = B;
        return;
    end
    allNames = unique([A.Properties.VariableNames, B.Properties.VariableNames], 'stable');
    A = add_missing_vars_(A, allNames);
    B = add_missing_vars_(B, allNames);
    out = [A(:, allNames); B(:, allNames)];
end

function T = add_missing_vars_(T, names)
    for i = 1:numel(names)
        n = names{i};
        if ~ismember(n, T.Properties.VariableNames)
            T.(n) = strings(height(T), 1);
        end
    end
end

function x = to_double_col_(x)
    if isnumeric(x)
        x = double(x(:));
    elseif islogical(x)
        x = double(x(:));
    else
        x = str2double(string(x(:)));
    end
end

function x = to_logical_col_(x)
    if islogical(x)
        x = x(:);
    elseif isnumeric(x)
        x = x(:) ~= 0;
    else
        s = lower(strtrim(string(x(:))));
        x = ismember(s, ["true", "1", "yes", "y", "pass", "passed"]);
    end
end
