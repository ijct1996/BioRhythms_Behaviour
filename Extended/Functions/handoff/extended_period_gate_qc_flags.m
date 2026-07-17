function Q = extended_period_gate_qc_flags(M)
%EXTENDED_PERIOD_GATE_QC_FLAGS Rows that failed CarryForward or need review flags.

    if isempty(M)
        Q = table();
        return;
    end
    flagMask = ~M.CarryForward | M.HarmonicSensitive12hFlag | ...
        contains(M.FullLadderSensitivityStatus, "sensitive", 'IgnoreCase', true);
    Q = M(flagMask, :);
end
