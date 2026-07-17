function Honly = extended_period_gate_hsub_only(Raw, HSub, tolLog2, primaryMode)
%EXTENDED_PERIOD_GATE_HSUB_ONLY Eligible primary HSub features without Raw match.

    hdr = {'File', 'SignalID', 'ConditionParsed', 'Photoperiod_h', 'Phase', 'BandName', ...
        'PrimaryHSubMode', 'HSubGlobalCandidateRow', 'HSubCandidateID', 'HSubPeriod_h', ...
        'HSubRidgeCoverageFrac', 'HSubCOIValidFrac', 'HSubPassQC', 'Reason'};
    rows = {};

    if isempty(HSub) || height(HSub) == 0
        Honly = cell2table(cell(0, numel(hdr)), 'VariableNames', hdr);
        return;
    end

    for i = 1:height(HSub)
        h = HSub(i, :);
        if ~logical(h.EligibleHSub), continue; end

        rAny = Raw(extended_period_gate_same_key(Raw, h), :);
        if isempty(rAny)
            reason = "No_raw_candidate";
            isHOnly = true;
        else
            if ismember('EligibleRaw', rAny.Properties.VariableNames)
                rAny = rAny(rAny.EligibleRaw, :);
            end
            if isempty(rAny)
                reason = "Raw_candidate_QC_fail";
                isHOnly = true;
            else
                d = abs(log2(rAny.MedianRidgePeriod_h ./ h.MedianRidgePeriod_h));
                isHOnly = all(d > tolLog2);
                reason = "No_raw_period_match_within_tolerance";
            end
        end

        if isHOnly
            rows(end+1, :) = {string(h.File), string(h.SignalID), string(h.ConditionParsed), ...
                h.Photoperiod_h, string(h.Phase), string(h.BandName), string(primaryMode), ...
                h.GlobalCandidateRow, string(h.CandidateID), h.MedianRidgePeriod_h, ...
                h.RidgeCoverageFrac, h.COIValidFrac, logical(h.PassQC), reason}; %#ok<AGROW>
        end
    end

    if isempty(rows)
        Honly = cell2table(cell(0, numel(hdr)), 'VariableNames', hdr);
    else
        Honly = cell2table(rows, 'VariableNames', hdr);
    end
end
