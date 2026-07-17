function M = extended_period_gate_add_fl_sensitivity(M, flModes)
%EXTENDED_PERIOD_GATE_ADD_FL_SENSITIVITY Full-Ladder sensitivity status on matched rows.

    status = strings(height(M), 1);
    for i = 1:height(M)
        if ~M.CarryForward(i)
            status(i) = "Not_applicable_not_carried_forward";
            continue;
        end
        anyExists = false;
        anyMatched = false;
        for j = 1:numel(flModes)
            prefix = matlab.lang.makeValidName(char(flModes(j)));
            existsCol = [prefix '_Exists'];
            matchedCol = [prefix '_Matched'];
            if ismember(existsCol, M.Properties.VariableNames)
                anyExists = anyExists || M.(existsCol)(i);
            end
            if ismember(matchedCol, M.Properties.VariableNames)
                anyMatched = anyMatched || M.(matchedCol)(i);
            end
        end
        if anyMatched
            status(i) = "Selective_supported_and_FullLadder_supported";
        elseif anyExists
            status(i) = "Selective_supported_FullLadder_sensitive";
        else
            status(i) = "Selective_supported_NoFullLadderCandidate";
        end
    end
    M.FullLadderSensitivityStatus = status;
end
