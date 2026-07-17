function [periodTable, ridgePhaseLong] = extended_ridge_build_tables( ...
    fileStem, signalID, sourceTag, hsubModeTag, photoperiod_h, ...
    bandNames, bandsMat, time_day, ZT_hr, lightStateStr, phaseMasks, bandTS, cfg)
%EXTENDED_RIDGE_BUILD_TABLES PeriodCandidates rows + RidgePhase_Long for one signal.

    minCov = cfg.ridgeHandoff.minRidgeCoverage;
    minCOI = cfg.ridgeHandoff.minCOIValidFrac;
    condParsed = extended_ridge_parse_condition(signalID);

    periodRows = {};
    phaseRows = {};
    nT = numel(time_day);
    phaseList = {'All'};
    if isfield(phaseMasks, 'Light')
        phaseList = {'All', 'Light', 'Dark'};
    end

    periodHdr = {'File', 'SignalID', 'ConditionParsed', 'Source', 'HSubResidualMode', ...
        'Photoperiod_h', 'Phase', 'BandName', 'CandidateID', 'CandidateRank', ...
        'MedianRidgePeriod_h', 'IQR_RidgePeriod_h', 'MeanBandPower_log10', 'SDBandPower_log10', ...
        'MeanRidgePower_log10', 'SDRidgePower_log10', 'RidgeCoverageFrac', 'COIValidFrac', ...
        'ValidPointCount', 'TotalPointCount', 'PassQC', 'QCReason'};

    phaseHdr = {'File', 'SignalID', 'ConditionParsed', 'Source', 'HSubResidualMode', ...
        'Photoperiod_h', 'BandName', 'CandidateID', 'Time_days', 'ZT_hours', ...
        'LightStateValue', 'Phase', 'RidgePeriod_h', 'RidgePower_log10', ...
        'RidgePhase_rad', 'ValidFlag'};

    for p = 1:numel(phaseList)
        phaseTag = phaseList{p};
        mask = extended_ridge_resize_mask(phaseMasks.(phaseTag), nT);
        totalN = sum(mask);

        for b = 1:numel(bandNames)
            bn = bandNames{b};
            vOK = logical(bandTS.ValidFlag(b, :).') & mask(:);
            rp = bandTS.RidgePeriod(b, :).';
            bp = bandTS.BandPower(b, :).';
            rpw = bandTS.RidgePower(b, :).';
            ph = bandTS.RidgePhase(b, :).';

            validN = sum(vOK & isfinite(rp));
            coverage = validN / max(totalN, 1);
            coiFrac = coverage;
            medRP = median(rp(vOK), 'omitnan');
            iqrRP = extended_ridge_iqr(rp(vOK));
            passQC = isfinite(medRP) && coverage >= minCov && coiFrac >= minCOI;
            if passQC
                qcReason = 'Pass';
            elseif ~isfinite(medRP)
                qcReason = 'NoFiniteRidgePeriod';
            elseif coverage < minCov
                qcReason = sprintf('LowRidgeCoverage_<%.3g', minCov);
            else
                qcReason = sprintf('LowCOIValidFrac_<%.3g', minCOI);
            end

            candID = extended_ridge_make_candidate_id(fileStem, signalID, sourceTag, ...
                hsubModeTag, phaseTag, bn, 1);
            periodRows(end+1, :) = {fileStem, signalID, condParsed, sourceTag, hsubModeTag, ...
                photoperiod_h, phaseTag, bn, candID, 1, medRP, iqrRP, ...
                mean(bp(vOK), 'omitnan'), std(bp(vOK), 0, 'omitnan'), ...
                mean(rpw(vOK), 'omitnan'), std(rpw(vOK), 0, 'omitnan'), ...
                coverage, coiFrac, validN, totalN, double(passQC), qcReason}; %#ok<AGROW>

            for t = 1:nT
                if ~mask(t), continue; end
                phaseRows(end+1, :) = {fileStem, signalID, condParsed, sourceTag, hsubModeTag, ...
                    photoperiod_h, bn, candID, time_day(t), ZT_hr(t), lightStateStr{t}, ...
                    phaseTag, rp(t), rpw(t), ph(t), double(bandTS.ValidFlag(b, t))}; %#ok<AGROW>
            end
        end
    end

    if isempty(periodRows)
        periodTable = cell2table(cell(0, numel(periodHdr)), 'VariableNames', periodHdr);
    else
        periodTable = cell2table(periodRows, 'VariableNames', periodHdr);
    end

    if isempty(phaseRows)
        ridgePhaseLong = cell2table(cell(0, numel(phaseHdr)), 'VariableNames', phaseHdr);
    else
        ridgePhaseLong = cell2table(phaseRows, 'VariableNames', phaseHdr);
        ridgePhaseLong = ridgePhaseLong(strcmpi(string(ridgePhaseLong.Phase), 'All'), :);
    end
end

function cond = extended_ridge_parse_condition(signalID)
    s = char(string(signalID));
    k = find(s == '_', 1, 'first');
    if isempty(k) || k == 1
        cond = s;
    else
        cond = s(1:k-1);
    end
end

function m2 = extended_ridge_resize_mask(m, N)
    m = logical(m(:));
    if numel(m) < N
        if isempty(m)
            m2 = false(N, 1);
        else
            m2 = [m; repmat(m(end), N - numel(m), 1)];
        end
    elseif numel(m) > N
        m2 = m(1:N);
    else
        m2 = m;
    end
end

function val = extended_ridge_iqr(x)
    x = x(isfinite(x));
    if isempty(x)
        val = NaN;
    else
        val = iqr(x);
    end
end
