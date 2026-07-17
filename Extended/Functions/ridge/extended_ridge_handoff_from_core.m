function summary = extended_ridge_handoff_from_core(entries, extHandoffDir, cfg)
%EXTENDED_RIDGE_HANDOFF_FROM_CORE Ridge-only WP handoff from Core Script 1–3 outputs.
%
%   Reads CoreSummary paths (raw xlsx, HSub residuals). Does NOT re-run HSub or
%   Core wavelet scalograms. Writes Kent-compatible WP_Summary__ / WP_TS__ mats.
%
%   Per-photoperiod failures (including LL edge cases) are skipped and logged;
%   remaining files still complete.

    if nargin < 3 || isempty(cfg)
        cfg = extended_defaults();
    end

    extended_period_gate_ensure_dir(extHandoffDir);
    logPath = fullfile(extHandoffDir, 'Logs', sprintf('Script4_RidgeHandoff_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    extended_period_gate_ensure_dir(fileparts(logPath));
    LOG = fopen(logPath, 'w');
    cleanupObj = onCleanup(@() extended_period_gate_fclose_if_open(LOG)); %#ok<NASGU>

    bandsMat = [reshape(cfg.bands.CR, 1, 2); cfg.bands.UR];
    bandNames = cellstr(cfg.bands.allNames);

    hsubModes = unique([string(cfg.hsub.primaryMode); ...
        string(cfg.hsub.secondaryModes(:)); ...
        string(cfg.hsub.sensitivityModes(:))], 'stable');

    loadLog = {};
    skipLog = {};

    summary = struct();
    summary.extHandoffDir = extHandoffDir;
    summary.logPath = logPath;

    for i = 1:numel(entries)
        entry = entries(i);
        fileStem = char(string(entry.fileStem));
        lightH = entry.lightHours;
        fprintf('Script 4 ridge handoff: %s (L%.0f h)\n', fileStem, lightH);
        extended_period_gate_log(LOG, 'Processing %s (L%.0f)', fileStem, lightH);

        try
            fileInfo = process_one_entry_(entry, fileStem, lightH, extHandoffDir, ...
                bandsMat, bandNames, hsubModes, cfg, LOG);
            loadLog{end+1} = fileInfo; %#ok<AGROW>
            fprintf('Wrote %s (%d candidates, %d phase rows)\n', ...
                fileStem, fileInfo.nCandidates, fileInfo.nPhaseRows);
            extended_period_gate_log(LOG, 'Wrote %s (%d candidates, %d phase rows)', ...
                fileStem, fileInfo.nCandidates, fileInfo.nPhaseRows);
        catch ME
            skipLog{end+1} = struct( ...
                'fileStem', fileStem, ...
                'lightHours', lightH, ...
                'message', ME.message, ...
                'identifier', ME.identifier); %#ok<AGROW>
            warning('extended_ridge_handoff_from_core:SkippedFile', ...
                'Skipping %s (L%.0f h): %s', fileStem, lightH, ME.message);
            extended_period_gate_log(LOG, 'SKIPPED %s (L%.0f): %s', fileStem, lightH, ME.message);
        end
    end

    if isempty(loadLog)
        summary.files = struct([]);
    else
        summary.files = [loadLog{:}];
    end

    if isempty(skipLog)
        summary.skipped = struct([]);
    else
        summary.skipped = [skipLog{:}];
    end

    fprintf('Script 4 ridge handoff complete: %s\n', extHandoffDir);
    fprintf('  Succeeded: %d | Skipped: %d\n', numel(loadLog), numel(skipLog));
    extended_period_gate_log(LOG, 'Complete. Succeeded=%d Skipped=%d', numel(loadLog), numel(skipLog));

    if isempty(loadLog)
        error('extended_ridge_handoff_from_core:NoFiles', ...
            'No WP handoff files were written. See log: %s', logPath);
    end
end

%% ------------------------------------------------------------------------
function fileInfo = process_one_entry_(entry, fileStem, lightH, extHandoffDir, ...
        bandsMat, bandNames, hsubModes, cfg, LOG)

    if ~isfield(entry, 'paths') || ~isfield(entry.paths, 'raw') || ~isfile(entry.paths.raw)
        error('extended_ridge_handoff_from_core:NoRaw', 'Missing raw path for %s', fileStem);
    end

    [rawTbl, rawMeta] = read_behav_excel(entry.paths.raw);
    time_hr = rawTbl{:, rawMeta.timeIdx};
    time_day = time_hr / 24;
    ZT_hr = mod(time_hr, 24);
    lightVec = rawTbl{:, rawMeta.lightIdx};
    [~, lightStateStr, phaseMasks] = extended_ridge_phase_masks(time_hr, lightVec);

    photoperiod_h = lightH;
    if ~isfinite(photoperiod_h)
        photoperiod_h = rawMeta.lightDurationHours;
    end

    % LL / constant light: still analysed; note in log (no true dark phase)
    if isfinite(photoperiod_h) && photoperiod_h >= cfg.ll.photoperiodValue
        extended_period_gate_log(LOG, ...
            'Note: %s is LL/constant-light (L>=%.0f). Ridge handoff proceeds; true LD/DL transitions deferred to Script 5.', ...
            fileStem, cfg.ll.photoperiodValue);
    end

    nSamples = height(rawTbl);
    FB = wavelet_make_filterbank(nSamples, cfg.samplingMinutes, core_defaults());

    periodAll = {};
    phaseAll = {};

    groups = entry.groups;
    if isempty(groups)
        error('extended_ridge_handoff_from_core:NoGroups', 'No groups for %s', fileStem);
    end

    for g = 1:numel(groups)
        if isfield(entry.meta, 'varNames') && ~isempty(entry.meta.varNames)
            grp = resolve_groups_by_names(entry.meta.varNames, groups(g));
        else
            grp = groups(g);
        end
        if isempty(grp.colIdx), continue; end

        grpName = grp.name;
        for m = 1:numel(grp.colIdx)
            colIdx = grp.colIdx(m);
            mouseName = rawMeta.varNames{colIdx};
            signalID = sprintf('%s_%s', grpName, mouseName);

            % --- RAW ridge ---
            signal = wavelet_prepare_signal(rawTbl, colIdx);
            [wt, periods_hours, coi_hours] = wavelet_compute_cwt(signal, FB);
            logPow = log10(max(abs(wt).^2, eps));
            bandTS = extended_ridge_compute_band_ts(periods_hours, logPow, bandsMat, coi_hours, wt);
            [pTab, phTab] = extended_ridge_build_tables(fileStem, signalID, 'Raw', 'NA', ...
                photoperiod_h, bandNames, bandsMat, time_day, ZT_hr, lightStateStr, phaseMasks, bandTS, cfg);
            periodAll{end+1} = pTab; %#ok<AGROW>
            phaseAll{end+1} = phTab; %#ok<AGROW>

            % --- HSub ridge modes ---
            if isfield(entry, 'paths') && isfield(entry.paths, 'hsub')
                for hm = 1:numel(hsubModes)
                    mode = hsubModes(hm);
                    sigH = extended_ridge_read_hsub_signal(entry.paths.hsub, fileStem, mode, mouseName, cfg);
                    if isempty(sigH), continue; end
                    sigH = extended_wavelet_prepare_vector(sigH);
                    [wtH, periods_hours, coi_hours] = wavelet_compute_cwt(sigH, FB);
                    logPowH = log10(max(abs(wtH).^2, eps));
                    bandTSH = extended_ridge_compute_band_ts(periods_hours, logPowH, bandsMat, coi_hours, wtH);
                    [pTabH, phTabH] = extended_ridge_build_tables(fileStem, signalID, 'HSub', char(mode), ...
                        photoperiod_h, bandNames, bandsMat, time_day, ZT_hr, lightStateStr, phaseMasks, bandTSH, cfg);
                    periodAll{end+1} = pTabH; %#ok<AGROW>
                    phaseAll{end+1} = phTabH; %#ok<AGROW>
                end
            end
        end
    end

    PeriodCandidates_Long = extended_ridge_vertcat_tables(periodAll);
    RidgePhase_Long = extended_ridge_vertcat_tables(phaseAll);

    if isempty(PeriodCandidates_Long) || height(PeriodCandidates_Long) == 0
        error('extended_ridge_handoff_from_core:NoCandidates', ...
            'No period candidates produced for %s', fileStem);
    end

    pkgS = struct();
    pkgS.meta = struct('script', 'extended_script4', 'fileStem', fileStem, ...
        'source', 'Core handoff ridge-only', 'timestamp', datetime('now'));
    pkgS.file = struct('FileStem', fileStem, 'InputFile', entry.paths.raw);
    pkgS.tables = struct();
    pkgS.tables.PeriodCandidates_Long = PeriodCandidates_Long;

    pkgTS = struct();
    pkgTS.meta = pkgS.meta;
    pkgTS.file = pkgS.file;
    pkgTS.tables = struct();
    pkgTS.tables.RidgePhase_Long = RidgePhase_Long;

    sumPath = fullfile(extHandoffDir, sprintf('WP_Summary__%s.mat', fileStem));
    tsPath = fullfile(extHandoffDir, sprintf('WP_TS__%s.mat', fileStem));
    save(sumPath, 'pkgS', '-v7.3');
    save(tsPath, 'pkgTS', '-v7.3');

    fileInfo = struct( ...
        'fileStem', fileStem, ...
        'lightHours', photoperiod_h, ...
        'summaryPath', sumPath, ...
        'tsPath', tsPath, ...
        'nCandidates', height(PeriodCandidates_Long), ...
        'nPhaseRows', height(RidgePhase_Long));
end

function sig = extended_ridge_read_hsub_signal(hsubRoot, fileStem, mode, mouseName, cfg)
    sig = [];
    mode = string(mode);
    if ~isfield(cfg.hsub, 'residualArms') || ~isfield(cfg.hsub.residualArms, char(mode))
        return;
    end
    armKey = cfg.hsub.residualArms.(char(mode));
    resType = cfg.hsub.residualTypes.(char(mode));
    resFile = fullfile(hsubRoot, fileStem, 'TimeSeries', 'Residual', ...
        sprintf('Residual_%s_%s_%s.xlsx', resType, armKey, fileStem));
    if ~isfile(resFile), return; end

    resTbl = readtable(resFile);
    resVar = matlab.lang.makeValidName(sprintf('Residual_%s_%s', resType, mouseName));
    if ~ismember(resVar, resTbl.Properties.VariableNames)
        resVar = matlab.lang.makeValidName(sprintf('Residual_Selective_%s', mouseName));
    end
    if ~ismember(resVar, resTbl.Properties.VariableNames)
        return;
    end
    sig = resTbl.(resVar);
end

function T = extended_ridge_vertcat_tables(cells)
    T = table();
    for k = 1:numel(cells)
        if isempty(cells{k}) || ~istable(cells{k}), continue; end
        if isempty(T)
            T = cells{k};
        else
            allVars = unique([T.Properties.VariableNames, cells{k}.Properties.VariableNames], 'stable');
            T = extended_ridge_add_missing(T, allVars);
            B = extended_ridge_add_missing(cells{k}, allVars);
            T = [T(:, allVars); B(:, allVars)]; %#ok<AGROW>
        end
    end
end

function T = extended_ridge_add_missing(T, allVars)
    for i = 1:numel(allVars)
        v = allVars{i};
        if ~ismember(v, T.Properties.VariableNames)
            T.(v) = repmat(missing, height(T), 1);
        end
    end
end
