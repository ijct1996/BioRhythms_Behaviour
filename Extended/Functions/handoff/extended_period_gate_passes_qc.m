function tf = extended_period_gate_passes_qc(T, minCoverage, minCOI, requirePassQC)
%EXTENDED_PERIOD_GATE_PASSES_QC Ridge / COI / PassQC eligibility for candidates.

    tf = true(height(T), 1);
    if requirePassQC && ismember('PassQC', T.Properties.VariableNames)
        tf = tf & logical(T.PassQC);
    end
    if ismember('RidgeCoverageFrac', T.Properties.VariableNames)
        tf = tf & T.RidgeCoverageFrac >= minCoverage;
    end
    if ismember('COIValidFrac', T.Properties.VariableNames)
        tf = tf & T.COIValidFrac >= minCOI;
    end
    if ismember('MedianRidgePeriod_h', T.Properties.VariableNames)
        tf = tf & isfinite(T.MedianRidgePeriod_h) & T.MedianRidgePeriod_h > 0;
    end
end
