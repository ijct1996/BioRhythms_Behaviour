function M = extended_period_gate_add_final_class(M)
%EXTENDED_PERIOD_GATE_ADD_FINAL_CLASS FinalValidationClass labels for matched rows.

    cls = string(M.MatchStatus);
    for i = 1:height(M)
        if M.CarryForward(i)
            cls(i) = "HSub_supported_primary";
            if strcmpi(M.FullLadderSensitivityStatus(i), "Selective_supported_FullLadder_sensitive")
                cls(i) = cls(i) + ";FullLadder_sensitive";
            end
            if M.HarmonicSensitive12hFlag(i)
                cls(i) = cls(i) + ";Harmonic_sensitive_12h_region";
            end
        else
            cls(i) = string(M.MatchStatus(i));
            if M.HarmonicSensitive12hFlag(i)
                cls(i) = cls(i) + ";Harmonic_sensitive_12h_region";
            end
        end
    end
    M.FinalValidationClass = cls;
end
