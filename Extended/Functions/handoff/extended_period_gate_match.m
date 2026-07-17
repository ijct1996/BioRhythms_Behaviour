function M = extended_period_gate_match(Raw, HSub, primaryMode, tolLog2, tolFrac, harmonicBands)
%EXTENDED_PERIOD_GATE_MATCH Primary Raw vs SEL_P360 (or configured) period matches.
%
%   CarryForward = true only when eligible Raw finds eligible HSub within tolLog2.
%   tolFrac is retained for API symmetry with legacy (log2 tolerance drives matching).

    %#ok<INUSD> tolFrac
    varNames = {'File', 'SignalID', 'ConditionParsed', 'Photoperiod_h', 'Phase', 'BandName', ...
        'RawGlobalCandidateRow', 'RawCandidateID', 'RawCandidateRank', 'RawPeriod_h', 'RawIQR_h', ...
        'RawMeanBandPower_log10', 'RawMeanRidgePower_log10', 'RawRidgeCoverageFrac', ...
        'RawCOIValidFrac', 'RawPassQC', ...
        'PrimaryHSubMode', 'PrimaryHSubGlobalCandidateRow', 'PrimaryHSubCandidateID', ...
        'PrimaryHSubPeriod_h', 'PrimaryHSubIQR_h', ...
        'PrimaryHSubMeanRidgePower_log10', 'PrimaryHSubRidgeCoverageFrac', ...
        'PrimaryHSubCOIValidFrac', 'PrimaryHSubPassQC', ...
        'AbsLog2PeriodDiff', 'PeriodDiffPercent', 'WithinTolerance', 'MatchStatus', ...
        'CarryForward', 'HarmonicSensitive12hFlag', 'QCFlag', 'SourcePackage'};

    rows = cell(height(Raw), numel(varNames));

    for i = 1:height(Raw)
        r = Raw(i, :);
        rawOK = logical(r.EligibleRaw);
        rawP = r.MedianRidgePeriod_h;
        harmonicFlag = ismember(string(r.BandName), harmonicBands);

        hsubAny = table();
        hsubElig = table();
        if ~isempty(HSub)
            keyMask = extended_period_gate_same_key(HSub, r);
            hsubAny = HSub(keyMask, :);
            hsubElig = HSub(keyMask & HSub.EligibleHSub, :);
        end

        best = table();
        absLog2Diff = NaN;
        pctDiff = NaN;
        withinTol = false;

        if rawOK && ~isempty(hsubElig)
            d = abs(log2(rawP ./ hsubElig.MedianRidgePeriod_h));
            [absLog2Diff, idx] = min(d);
            best = hsubElig(idx, :);
            pctDiff = 100 * abs(rawP - best.MedianRidgePeriod_h) / rawP;
            withinTol = absLog2Diff <= tolLog2;
        elseif ~rawOK && ~isempty(hsubAny)
            d = abs(log2(rawP ./ hsubAny.MedianRidgePeriod_h));
            [absLog2Diff, idx] = min(d);
            best = hsubAny(idx, :);
            pctDiff = 100 * abs(rawP - best.MedianRidgePeriod_h) / rawP;
            withinTol = false;
        elseif rawOK && isempty(hsubElig) && ~isempty(hsubAny)
            d = abs(log2(rawP ./ hsubAny.MedianRidgePeriod_h));
            [absLog2Diff, idx] = min(d);
            best = hsubAny(idx, :);
            pctDiff = 100 * abs(rawP - best.MedianRidgePeriod_h) / rawP;
            withinTol = false;
        end

        if ~rawOK
            status = "Fail_raw_QC";
            carry = false;
            qc = "Raw_QC_fail";
        elseif isempty(hsubAny)
            status = "Raw_only_not_carried_forward";
            carry = false;
            qc = "No_primary_HSub_candidate";
        elseif isempty(hsubElig)
            status = "Raw_present_primary_HSub_QC_fail";
            carry = false;
            qc = "Primary_HSub_QC_fail";
        elseif withinTol
            status = "HSub_supported_primary";
            carry = true;
            qc = "Pass";
        else
            status = "Raw_present_HSub_period_mismatch";
            carry = false;
            qc = "Period_mismatch";
        end

        if harmonicFlag
            qc = qc + ";Harmonic_sensitive_12h_region";
        end

        [hID, hRow, hP, hIQR, hPow, hCov, hCOI, hPass] = empty_hsub_vals_();
        if ~isempty(best)
            hID   = string(best.CandidateID);
            hRow  = best.GlobalCandidateRow;
            hP    = best.MedianRidgePeriod_h;
            hIQR  = best.IQR_RidgePeriod_h;
            hPow  = best.MeanRidgePower_log10;
            hCov  = best.RidgeCoverageFrac;
            hCOI  = best.COIValidFrac;
            hPass = logical(best.PassQC);
        end

        rows(i, :) = {string(r.File), string(r.SignalID), string(r.ConditionParsed), ...
            r.Photoperiod_h, string(r.Phase), string(r.BandName), ...
            r.GlobalCandidateRow, string(r.CandidateID), r.CandidateRank, ...
            r.MedianRidgePeriod_h, r.IQR_RidgePeriod_h, ...
            r.MeanBandPower_log10, r.MeanRidgePower_log10, r.RidgeCoverageFrac, ...
            r.COIValidFrac, rawOK, ...
            string(primaryMode), hRow, hID, hP, hIQR, hPow, hCov, hCOI, hPass, ...
            absLog2Diff, pctDiff, withinTol, status, carry, harmonicFlag, qc, ...
            string(r.SourcePackage)};
    end

    M = cell2table(rows, 'VariableNames', varNames);
end

function [hID, hRow, hP, hIQR, hPow, hCov, hCOI, hPass] = empty_hsub_vals_()
    hID = string(missing);
    hRow = NaN;
    hP = NaN;
    hIQR = NaN;
    hPow = NaN;
    hCov = NaN;
    hCOI = NaN;
    hPass = false;
end
