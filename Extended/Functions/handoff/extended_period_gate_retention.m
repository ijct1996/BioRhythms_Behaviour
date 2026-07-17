function S = extended_period_gate_retention(M, keyVars)
%EXTENDED_PERIOD_GATE_RETENTION Retention summary by Band / Photoperiod / both.

    keyVars = cellstr(keyVars);
    if isempty(M)
        S = table();
        return;
    end
    [U, ~, ic] = unique(M(:, keyVars), 'rows', 'stable');
    n = height(U);
    N_RawCandidates = zeros(n, 1);
    N_RawEligible = zeros(n, 1);
    N_CarryForward = zeros(n, 1);
    N_HarmonicSensitive12h = zeros(n, 1);
    RetentionPct_EligibleRaw = nan(n, 1);
    RetentionPct_AllRaw = nan(n, 1);

    for i = 1:n
        idx = ic == i;
        N_RawCandidates(i) = sum(idx);
        N_RawEligible(i) = sum(M.RawPassQC(idx));
        N_CarryForward(i) = sum(M.CarryForward(idx));
        N_HarmonicSensitive12h(i) = sum(M.HarmonicSensitive12hFlag(idx));
        if N_RawEligible(i) > 0
            RetentionPct_EligibleRaw(i) = 100 * N_CarryForward(i) / N_RawEligible(i);
        end
        if N_RawCandidates(i) > 0
            RetentionPct_AllRaw(i) = 100 * N_CarryForward(i) / N_RawCandidates(i);
        end
    end

    S = U;
    S.N_RawCandidates = N_RawCandidates;
    S.N_RawEligible = N_RawEligible;
    S.N_CarryForward = N_CarryForward;
    S.RetentionPct_EligibleRaw = RetentionPct_EligibleRaw;
    S.RetentionPct_AllRaw = RetentionPct_AllRaw;
    S.N_HarmonicSensitive12h = N_HarmonicSensitive12h;
end
