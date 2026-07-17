function M = extended_period_gate_add_mode_columns(M, Hm, mode, tolLog2)
%EXTENDED_PERIOD_GATE_ADD_MODE_COLUMNS Secondary / Full-Ladder mode match columns.

    prefix = matlab.lang.makeValidName(char(mode));
    existsCol = [prefix '_Exists'];
    matchedCol = [prefix '_Matched'];
    periodCol = [prefix '_Period_h'];
    diffCol = [prefix '_DiffPercent'];
    idCol = [prefix '_CandidateID'];
    eligibleCol = [prefix '_EligibleHSub'];

    M.(existsCol) = false(height(M), 1);
    M.(matchedCol) = false(height(M), 1);
    M.(periodCol) = nan(height(M), 1);
    M.(diffCol) = nan(height(M), 1);
    M.(idCol) = strings(height(M), 1);
    M.(eligibleCol) = false(height(M), 1);

    if isempty(Hm) || height(Hm) == 0
        return;
    end

    for i = 1:height(M)
        rawProxy = table();
        rawProxy.File = M.File(i);
        rawProxy.SignalID = M.SignalID(i);
        rawProxy.ConditionParsed = M.ConditionParsed(i);
        rawProxy.Photoperiod_h = M.Photoperiod_h(i);
        rawProxy.Phase = M.Phase(i);
        rawProxy.BandName = M.BandName(i);

        keyMask = extended_period_gate_same_key(Hm, rawProxy);
        hAny = Hm(keyMask, :);
        if ~isempty(hAny)
            M.(existsCol)(i) = true;
        else
            continue;
        end

        hElig = hAny;
        if ismember('EligibleHSub', hElig.Properties.VariableNames)
            hElig = hElig(hElig.EligibleHSub, :);
        end
        if isempty(hElig)
            continue;
        end

        M.(eligibleCol)(i) = true;
        rawP = M.RawPeriod_h(i);
        d = abs(log2(rawP ./ hElig.MedianRidgePeriod_h));
        [bestD, idx] = min(d);
        best = hElig(idx, :);
        M.(periodCol)(i) = best.MedianRidgePeriod_h;
        M.(diffCol)(i) = 100 * abs(rawP - best.MedianRidgePeriod_h) / rawP;
        M.(idCol)(i) = string(best.CandidateID);
        M.(matchedCol)(i) = bestD <= tolLog2;
    end
end
