function T = extended_period_gate_fl_sensitivity_table(M, flModes)
%EXTENDED_PERIOD_GATE_FL_SENSITIVITY_TABLE Slim Full-Ladder sensitivity view of Matched.

    baseVars = {'File', 'SignalID', 'ConditionParsed', 'Photoperiod_h', 'Phase', 'BandName', ...
        'RawCandidateID', 'RawPeriod_h', 'CarryForward', 'FullLadderSensitivityStatus', ...
        'FinalValidationClass'};
    extraVars = {};
    for j = 1:numel(flModes)
        prefix = matlab.lang.makeValidName(char(flModes(j)));
        extraVars = [extraVars, {[prefix '_Exists'], [prefix '_Matched'], ...
            [prefix '_Period_h'], [prefix '_DiffPercent'], [prefix '_CandidateID']}]; %#ok<AGROW>
    end
    vars = [baseVars, extraVars];
    vars = vars(ismember(vars, M.Properties.VariableNames));
    T = M(:, vars);
end
