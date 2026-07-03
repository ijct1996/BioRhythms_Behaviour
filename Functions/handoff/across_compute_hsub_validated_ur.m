function [bandTable, periodTable, validationLog] = across_compute_hsub_validated_ur(entries, cfg)
%ACROSS_COMPUTE_HSUB_VALIDATED_UR CWT on SEL_P360 residual; UR peaks/bands with harmonic filter.
%
%   Confirmatory Script 3 layer. Script 2 unchanged (RAW handoff preserved).
%   - CR_20_28 band power: RAW activity (cfg.hsubValidation.crSource)
%   - UR band power + UR period peaks: HSub residual (AnchorOK mice only)
%   - Peaks harmonically related to anchor P0 are excluded

    if nargin < 2, cfg = core_defaults(); end

    bandRows = {};
    periodRows = {};
    logRows = {};

    armKey = cfg.hsub.scalogramArm;

    for i = 1:numel(entries)
        entryBand = {};
        entryPeriod = {};
        entryLog = {};

        if ~isfield(entries(i), 'hsub') || ~entries(i).hsub.available
            warning('across_compute_hsub_validated_ur:NoHSub', ...
                'HSub summary missing for %s; skipped.', entries(i).fileStem);
            continue;
        end
        if ~isfield(entries(i), 'paths') || ~isfield(entries(i).paths, 'hsub')
            warning('across_compute_hsub_validated_ur:NoHSubPath', ...
                'paths.hsub missing for %s; skipped.', entries(i).fileStem);
            continue;
        end

        hsubRun = fullfile(entries(i).paths.hsub, entries(i).fileStem);
        resFile = fullfile(hsubRun, 'TimeSeries', 'Residual', ...
            sprintf('Residual_Selective_%s_%s.xlsx', armKey, entries(i).fileStem));
        if ~isfile(resFile)
            warning('across_compute_hsub_validated_ur:NoResidual', ...
                'Residual file missing for %s (re-run Script 1).', entries(i).fileStem);
            continue;
        end

        rawTbl = [];
        rawMeta = [];
        if strcmpi(cfg.hsubValidation.crSource, 'raw')
            if ~isfile(entries(i).paths.raw)
                error('across_compute_hsub_validated_ur:NoRaw', ...
                    'RAW file not found: %s', entries(i).paths.raw);
            end
            [rawTbl, rawMeta] = read_behav_excel(entries(i).paths.raw);
        end

        resTbl = readtable(resFile);
        anchorMap = across_hsub_build_anchor_map(entries(i).hsub.anchor);
        nSamples = height(resTbl);
        FB = wavelet_make_filterbank(nSamples, cfg.samplingMinutes, cfg);

        groups = entries(i).groups;
        if isempty(groups)
            warning('across_compute_hsub_validated_ur:NoGroups', ...
                'No groups in handoff for %s; skipped.', entries(i).fileStem);
            continue;
        end

        for g = 1:numel(groups)
            if isfield(entries(i).meta, 'varNames') && ~isempty(entries(i).meta.varNames)
                grp = resolve_groups_by_names(entries(i).meta.varNames, groups(g));
            else
                grp = groups(g);
            end
            if isempty(grp.colIdx)
                continue;
            end
            grpName = grp.name;
            colNames = grp.colNames;
            if isstring(colNames), colNames = cellstr(colNames); end
            if isempty(colNames) && ~isempty(grp.colIdx)
                colNames = entries(i).meta.varNames(grp.colIdx);
            end

            for m = 1:numel(colNames)
                colName = char(string(colNames{m}));

                anchor = across_hsub_lookup_anchor(anchorMap, colName);
                if ~anchor.ok
                    continue;
                end

                resVar = matlab.lang.makeValidName(['Residual_Selective_' colName]);
                if ~ismember(resVar, resTbl.Properties.VariableNames)
                    continue;
                end

                residual = resTbl.(resVar);
                residual = wavelet_prepare_signal_from_vector(residual);
                [wt, periods_hours, ~] = wavelet_compute_cwt(residual, FB);
                powerSpec = mean(abs(wt).^2, 2);

                bpStruct = wavelet_band_power(powerSpec, periods_hours, cfg);

                crPower = NaN;
                if ~isempty(rawTbl)
                    rawColIdx = find(strcmp(rawMeta.varNames, colName), 1);
                    if ~isempty(rawColIdx)
                        rawSig = wavelet_prepare_signal(rawTbl, rawColIdx);
                        [wtRaw, perRaw, ~] = wavelet_compute_cwt(rawSig, FB);
                        rawBands = wavelet_band_power(mean(abs(wtRaw).^2, 2), perRaw, cfg);
                        crPower = rawBands.CR_20_28;
                    end
                end

                signalLabel = sprintf('%s_%s', grpName, colName);
                bpRow = struct2table(bpStruct);
                bpRow.SignalID = string(signalLabel);
                bpRow.Group = string(grpName);
                bpRow.MouseColumn = string(colName);
                bpRow.LightDuration_h = entries(i).lightHours;
                bpRow.CR_20_28 = crPower;
                bpRow.AnchorOK = true(height(bpRow), 1);
                bpRow.P0_h = anchor.P0_h;
                bpRow.HSubValidated = true(height(bpRow), 1);
                bpRow.AnalysisNote = repmat("UR bands from SEL_P360 residual; CR from RAW", height(bpRow), 1);
                entryBand{end + 1} = bpRow; %#ok<AGROW>

                [peaks, peakLog] = across_hsub_find_validated_ur_peaks( ...
                    powerSpec, periods_hours, anchor.P0_h, cfg);
                if height(peaks) > 0
                    peaks.SignalID = repmat(string(signalLabel), height(peaks), 1);
                    peaks.Group = repmat(string(grpName), height(peaks), 1);
                    peaks.MouseColumn = repmat(string(colName), height(peaks), 1);
                    peaks.LightDuration_h = repmat(entries(i).lightHours, height(peaks), 1);
                    peaks.P0_h = repmat(anchor.P0_h, height(peaks), 1);
                    peaks.HSubValidated = true(height(peaks), 1);
                    entryPeriod{end + 1} = peaks; %#ok<AGROW>
                end
                if ~isempty(peakLog)
                    peakLog.SignalID = repmat(string(signalLabel), height(peakLog), 1);
                    peakLog.Group = repmat(string(grpName), height(peakLog), 1);
                    peakLog.MouseColumn = repmat(string(colName), height(peakLog), 1);
                    peakLog.LightDuration_h = repmat(entries(i).lightHours, height(peakLog), 1);
                    peakLog.P0_h = repmat(anchor.P0_h, height(peakLog), 1);
                    entryLog{end + 1} = peakLog; %#ok<AGROW>
                end
            end
        end

        bandRows = [bandRows, entryBand]; %#ok<AGROW>
        periodRows = [periodRows, entryPeriod]; %#ok<AGROW>
        logRows = [logRows, entryLog]; %#ok<AGROW>
    end

    bandTable = wavelet_vertcat_tables(bandRows);
    periodTable = wavelet_vertcat_tables(periodRows);
    validationLog = wavelet_vertcat_tables(logRows);

    if isempty(bandTable) || ~istable(bandTable)
        bandTable = empty_band_table();
    end
    if isempty(periodTable) || ~istable(periodTable)
        periodTable = empty_period_table();
    end
    if isempty(validationLog) || ~istable(validationLog)
        validationLog = empty_validation_log_table();
    end
end

function signal = wavelet_prepare_signal_from_vector(x)
    signal = x(:);
    if ~isnumeric(signal)
        signal = str2double(string(signal));
    end
    signal(~isfinite(signal)) = 0;
end

function anchorMap = across_hsub_build_anchor_map(anchorTable)
    anchorMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if isempty(anchorTable) || ~istable(anchorTable)
        return;
    end
    colVar = pick_var(anchorTable, {'Column', 'colName', 'Mouse'});
    okVar = pick_var(anchorTable, {'AnchorOK', 'anchorOK'});
    p0Var = pick_var(anchorTable, {'Period_hours', 'P0_h', 'Period_h'});
    for r = 1:height(anchorTable)
        key = char(string(anchorTable.(colVar)(r)));
        okVal = anchorTable.(okVar)(r);
        if islogical(okVal)
            ok = okVal;
        elseif isnumeric(okVal)
            ok = okVal ~= 0;
        else
            ok = any(strcmpi(string(okVal), ["true", "1", "yes"]));
        end
        p0 = anchorTable.(p0Var)(r);
        if isnumeric(p0), p0 = double(p0); else, p0 = str2double(string(p0)); end
        anchorMap(key) = struct('ok', ok, 'P0_h', p0);
    end
end

function anchor = across_hsub_lookup_anchor(anchorMap, colName)
    anchor = struct('ok', false, 'P0_h', NaN);
    if ~isKey(anchorMap, colName)
        return;
    end
    anchor = anchorMap(colName);
    anchor.ok = anchor.ok && isfinite(anchor.P0_h) && anchor.P0_h > 0;
end

function v = pick_var(tbl, candidates)
    for i = 1:numel(candidates)
        if ismember(candidates{i}, tbl.Properties.VariableNames)
            v = candidates{i};
            return;
        end
    end
    error('across_hsub_build_anchor_map:MissingColumn', ...
        'Expected one of: %s', strjoin(candidates, ', '));
end

function [peaks, logTable] = across_hsub_find_validated_ur_peaks(powerSpec, periods_hours, P0_h, cfg)
    urLim = cfg.hsubValidation.urPeriodRangeHours;
    urMask = periods_hours >= urLim(1) & periods_hours <= urLim(2);
    logP = log10(powerSpec + eps);
    logP_ur = logP;
    logP_ur(~urMask) = -inf;

    [pks, locs, w, prom] = findpeaks(logP_ur, 'MinPeakProminence', 0);
    peaks = table();
    logTable = table();

    if isempty(pks)
        return;
    end

    candPeriod = periods_hours(locs);
    isHarm = arrayfun(@(p) across_hsub_is_harmonic_period(p, P0_h, cfg), candPeriod);
    logTable = table(candPeriod, pks, prom, w, isHarm, ...
        'VariableNames', {'CandidatePeriod_h', 'PeakValue_log10', 'PeakProminence', ...
        'PeakWidth', 'IsHarmonicOfP0'});

    keep = ~isHarm;
    if ~any(keep)
        return;
    end

    pks = pks(keep);
    locs = locs(keep);
    prom = prom(keep);
    w = w(keep);
    candPeriod = candPeriod(keep);

    [~, ord] = sort(pks, 'descend');
    topN = cfg.hsubValidation.topNPeaks;
    ord = ord(1:min(topN, numel(ord)));

    peaks = table((1:numel(ord))', candPeriod(ord), pks(ord), prom(ord), w(ord), ...
        'VariableNames', {'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', ...
        'PeakProminence', 'PeakWidth'});
end

function isHarm = across_hsub_is_harmonic_period(period_h, P0_h, cfg)
    isHarm = false;
    if ~isfinite(period_h) || ~isfinite(P0_h) || P0_h <= 0
        return;
    end
    tol = cfg.hsubValidation.harmonicToleranceFrac;
    urLim = cfg.hsubValidation.urPeriodRangeHours;
    for k = 1:cfg.hsubValidation.maxHarmonicOrder
        target = P0_h / k;
        if target < urLim(1) || target > urLim(2)
            continue;
        end
        if abs(period_h - target) / target <= tol
            isHarm = true;
            return;
        end
    end
end

function T = empty_band_table()
    T = table(string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
        double.empty(0,1), double.empty(0,1), logical.empty(0,1), double.empty(0,1), ...
        logical.empty(0,1), string.empty(0,1), ...
        'VariableNames', {'SignalID', 'Group', 'MouseColumn', 'LightDuration_h', ...
        'CR_20_28', 'AnchorOK', 'P0_h', 'HSubValidated', 'AnalysisNote'});
end

function T = empty_period_table()
    T = table(double.empty(0,1), double.empty(0,1), double.empty(0,1), ...
        double.empty(0,1), double.empty(0,1), string.empty(0,1), string.empty(0,1), ...
        string.empty(0,1), double.empty(0,1), double.empty(0,1), logical.empty(0,1), ...
        'VariableNames', {'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', ...
        'PeakProminence', 'PeakWidth', 'SignalID', 'Group', 'MouseColumn', ...
        'LightDuration_h', 'P0_h', 'HSubValidated'});
end

function T = empty_validation_log_table()
    T = table(double.empty(0,1), double.empty(0,1), double.empty(0,1), ...
        double.empty(0,1), logical.empty(0,1), string.empty(0,1), string.empty(0,1), ...
        string.empty(0,1), double.empty(0,1), double.empty(0,1), ...
        'VariableNames', {'CandidatePeriod_h', 'PeakValue_log10', 'PeakProminence', ...
        'PeakWidth', 'IsHarmonicOfP0', 'SignalID', 'Group', 'MouseColumn', ...
        'LightDuration_h', 'P0_h'});
end
