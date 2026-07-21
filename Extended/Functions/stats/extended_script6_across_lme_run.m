function out = extended_script6_across_lme_run(extHandoffDir, cfg)
%EXTENDED_SCRIPT6_ACROSS_LME_RUN Local across-photoperiod LME / FDR (no Kent).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    handoffDir = resolve_handoff_dir_(extHandoffDir);
    outRoot = fullfile(handoffDir, 'AcrossPhotoperiod_LME');
    dirs = struct('Run', outRoot, 'Tables', fullfile(outRoot, 'Tables'), ...
        'Logs', fullfile(outRoot, 'Logs'));
    extended_period_gate_ensure_dir(dirs.Tables);
    extended_period_gate_ensure_dir(dirs.Logs);

    logPath = fullfile(dirs.Logs, sprintf('Script6_AcrossLME_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    LOG = fopen(logPath, 'w');
    cleanupLog = onCleanup(@() extended_period_gate_fclose_if_open(LOG)); %#ok<NASGU>

    fprintf('\n=== Extended Script 6: Across-photoperiod LME / FDR ===\n');
    fprintf('Handoff: %s\n', handoffDir);
    fprintf('Output:  %s\n', outRoot);
    log_line_(LOG, 'Handoff: %s', handoffDir);

    mapPath = fullfile(handoffDir, 'RawVsSelectiveHSub_PeriodValidation', 'HSubSupported_PeriodMap.mat');
    if ~isfile(mapPath)
        error('extended_script6_across_lme_run:MissingMap', ...
            'HSubSupported_PeriodMap.mat not found at %s. Run Script 4 first.', mapPath);
    end
    VM = load(mapPath, 'validationMap');
    CarryForward = VM.validationMap.CarryForward_Periods;
    fprintf('Carry-forward validated Raw UR candidates: %d\n', height(CarryForward));
    log_line_(LOG, 'CarryForward rows: %d', height(CarryForward));

    sumPaths = dir(fullfile(handoffDir, 'WP_Summary__*.mat'));
    if isempty(sumPaths)
        error('extended_script6_across_lme_run:NoSummaries', ...
            'No WP_Summary__*.mat under %s. Run Script 4 first.', handoffDir);
    end

    BCS = table();
    for i = 1:numel(sumPaths)
        p = fullfile(sumPaths(i).folder, sumPaths(i).name);
        S = load(p, 'pkgS');
        if ~isfield(S, 'pkgS') || isempty(S.pkgS)
            continue;
        end
        pkgS = S.pkgS;
        if isfield(pkgS, 'tables') && isfield(pkgS.tables, 'BandConditionSummary') ...
                && ~isempty(pkgS.tables.BandConditionSummary)
            T = pkgS.tables.BandConditionSummary;
        elseif isfield(pkgS, 'tables') && isfield(pkgS.tables, 'PeriodCandidates_Long') ...
                && ~isempty(pkgS.tables.PeriodCandidates_Long)
            T = extended_build_band_condition_summary(pkgS.tables.PeriodCandidates_Long);
            log_line_(LOG, 'Built BandConditionSummary from PeriodCandidates: %s (%d rows)', p, height(T));
        else
            log_line_(LOG, 'Skipped summary (no BCS or PeriodCandidates): %s', p);
            continue;
        end

        T.PackageID = repmat(i, height(T), 1);
        T.SummaryMat = repmat(string(p), height(T), 1);
        if isfield(pkgS, 'file') && isfield(pkgS.file, 'FileStem')
            T.FileStemPkg = repmat(string(pkgS.file.FileStem), height(T), 1);
        else
            T.FileStemPkg = repmat(missing, height(T), 1);
        end
        BCS = [BCS; T]; %#ok<AGROW>
    end

    if isempty(BCS)
        error('extended_script6_across_lme_run:NoBCS', ...
            'No BandConditionSummary rows could be built from WP_Summary packages.');
    end

    CR_BAND = string(cfg.bands.allNames(1));
    UR_BANDS = string(cfg.bands.UR_names);
    SRC_CR = "Raw";
    SRC_UR = "Raw";

    BCS_AllRows = annotate_bcs_with_validation_(BCS, CarryForward, UR_BANDS, SRC_UR, true);
    BCS_used = filter_bcs_for_analysis_(BCS_AllRows, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, true);
    BCS_used.Photoperiod_h = double(BCS_used.Photoperiod_h);
    fprintf('BandConditionSummary rows: %d | analysis-used: %d\n', height(BCS_AllRows), height(BCS_used));
    log_line_(LOG, 'BCS all=%d used=%d', height(BCS_AllRows), height(BCS_used));

    if isempty(BCS_used)
        error('extended_script6_across_lme_run:EmptyBCS', ...
            'No rows after validated-Raw filtering. Check CarryForward_Periods keys.');
    end

    xlsxOut = fullfile(dirs.Tables, 'AcrossPhotoperiod_Outputs.xlsx');
    if isfile(xlsxOut), delete(xlsxOut); end
    writetable(BCS_AllRows, xlsxOut, 'Sheet', 'BandConditionSummary_All');
    writetable(BCS_used, xlsxOut, 'Sheet', 'BandCondSummary_AnalysisUsed');
    writetable(CarryForward, xlsxOut, 'Sheet', 'ValidationMap_Used');

    [Tpair, TpairSummary] = build_cr_ur_pairs_(BCS_used, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
    Tabs = build_absolute_summary_(BCS_used, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
    writetable(Tpair, xlsxOut, 'Sheet', 'CR_UR_Pairs_PerMouse');
    writetable(TpairSummary, xlsxOut, 'Sheet', 'CR_UR_Pairs_Summary');
    writetable(Tabs, xlsxOut, 'Sheet', 'AbsolutePower_Summary');

    lmeResult = run_lme_block_(BCS_used, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, cfg, LOG);
    if lmeResult.success
        lmeXlsx = fullfile(dirs.Tables, 'AcrossPhotoperiod_LME_Outputs.xlsx');
        if isfile(lmeXlsx), delete(lmeXlsx); end
        writetable(lmeResult.inferenceRaw, lmeXlsx, 'Sheet', 'LME_Inference_Raw');
        writetable(lmeResult.inferenceFdr, lmeXlsx, 'Sheet', 'LME_Inference_BH_FDR');
        fprintf('LME workbook: %s\n', lmeXlsx);
        log_line_(LOG, 'LME workbook: %s', lmeXlsx);
    else
        warning('extended_script6_across_lme_run:LmeSkipped', ...
            'No LME models converged. Summary tables were still written to %s', xlsxOut);
    end

    out = struct();
    out.handoffDir = handoffDir;
    out.outRoot = outRoot;
    out.xlsxOut = xlsxOut;
    out.nBCS = height(BCS_used);
    out.lme = lmeResult;
    out.logPath = logPath;

    fprintf('\nExtended Script 6 complete.\n');
    fprintf('  Tables: %s\n', dirs.Tables);
end

%% ------------------------------------------------------------------------
function handoffDir = resolve_handoff_dir_(extHandoffDir)
    if nargin >= 1 && ~isempty(extHandoffDir)
        handoffDir = char(string(extHandoffDir));
        if ~isfolder(handoffDir)
            error('extended_script6_across_lme_run:BadHandoff', 'Not a folder: %s', handoffDir);
        end
        return;
    end
    rootSel = uigetdir(pwd, 'Select AcrossPhotoperiod_Input (folder with WP_Summary__*.mat)');
    if isequal(rootSel, 0)
        error('extended_script6_across_lme_run:NoFolder', 'No handoff folder selected.');
    end
    handoffDir = char(rootSel);
end

function log_line_(LOG, fmt, varargin)
    if isempty(LOG) || LOG < 0
        return;
    end
    fprintf(LOG, ['[' datestr(now, 'yyyy-mm-dd HH:MM:SS') '] ' fmt '\n'], varargin{:});
end

function [Tpair, Tsum] = build_cr_ur_pairs_(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    CR = BCS(BCS.Source == SRC_CR & BCS.BandName == CR_BAND, :);
    CR = CR(:, {'File', 'SignalID', 'Photoperiod_h', 'Phase', 'MeanBandPower_log10', 'FileStemPkg'});
    CR.Properties.VariableNames{'MeanBandPower_log10'} = 'CR_Log10';

    UR = BCS(BCS.Source == SRC_UR & ismember(BCS.BandName, UR_BANDS), :);
    UR = UR(:, {'File', 'SignalID', 'Photoperiod_h', 'Phase', 'BandName', 'MeanBandPower_log10', 'FileStemPkg'});
    UR.Properties.VariableNames{'BandName'} = 'UR_Band';
    UR.Properties.VariableNames{'MeanBandPower_log10'} = 'UR_Log10';

    keys = {'File', 'SignalID', 'Photoperiod_h', 'Phase', 'FileStemPkg'};
    Tpair = innerjoin(UR, CR, 'Keys', keys);
    Tpair.Delta_log10 = Tpair.UR_Log10 - Tpair.CR_Log10;
    Tpair.Ratio = 10 .^ Tpair.Delta_log10;

    G = findgroups(Tpair.Photoperiod_h, Tpair.Phase, Tpair.UR_Band);
    meanCR = splitapply(@(x) mean(x, 'omitnan'), Tpair.CR_Log10, G);
    sdCR = splitapply(@(x) std(x, 'omitnan'), Tpair.CR_Log10, G);
    meanUR = splitapply(@(x) mean(x, 'omitnan'), Tpair.UR_Log10, G);
    sdUR = splitapply(@(x) std(x, 'omitnan'), Tpair.UR_Log10, G);
    meanD = splitapply(@(x) mean(x, 'omitnan'), Tpair.Delta_log10, G);
    sdD = splitapply(@(x) std(x, 'omitnan'), Tpair.Delta_log10, G);
    [pp, ph, ub] = splitapply(@(a, b, c) deal(a(1), b(1), c(1)), ...
        Tpair.Photoperiod_h, Tpair.Phase, Tpair.UR_Band, G);

    Tsum = table(pp, ph, ub, meanCR, sdCR, meanUR, sdUR, meanD, sdD, ...
        'VariableNames', {'Photoperiod_h', 'Phase', 'UR_Band', ...
        'Mean_CR_Log10', 'SD_CR_Log10', 'Mean_UR_Log10', 'SD_UR_Log10', ...
        'Mean_Delta_log10', 'SD_Delta_log10'});
    Tsum = sortrows(Tsum, {'Phase', 'UR_Band', 'Photoperiod_h'});
end

function Tabs = build_absolute_summary_(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    bandsAll = [CR_BAND UR_BANDS];
    rows = cell(0, 5);
    for b = 1:numel(bandsAll)
        bn = bandsAll(b);
        if bn == CR_BAND
            D = BCS(BCS.Source == SRC_CR & BCS.BandName == bn, :);
        else
            D = BCS(BCS.Source == SRC_UR & BCS.BandName == bn, :);
        end
        if isempty(D), continue; end
        G = findgroups(D.Photoperiod_h, D.Phase);
        pp = splitapply(@(x) x(1), D.Photoperiod_h, G);
        ph = splitapply(@(x) x(1), D.Phase, G);
        m = splitapply(@(x) mean(x, 'omitnan'), D.MeanBandPower_log10, G);
        s = splitapply(@(x) std(x, 'omitnan'), D.MeanBandPower_log10, G);
        n = numel(m);
        add = [num2cell(pp(:)) cellstr(ph(:)) cellstr(repmat(string(bn), n, 1)) num2cell(m(:)) num2cell(s(:))];
        rows = [rows; add]; %#ok<AGROW>
    end
    Tabs = cell2table(rows, 'VariableNames', {'Photoperiod_h', 'Phase', 'BandName', 'Mean_Log10', 'SD_Log10'});
    Tabs.Phase = string(Tabs.Phase);
    Tabs.BandName = string(Tabs.BandName);
    Tabs = sortrows(Tabs, {'Phase', 'BandName', 'Photoperiod_h'});
end

function TBP = build_bandpower_table_(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    rows = cell(0, 7);
    hdr = {'File', 'SignalID', 'Photoperiod_h', 'Phase', 'BandName', 'Source', 'MeanBandPower_log10'};
    CR = BCS(BCS.Source == SRC_CR & BCS.BandName == CR_BAND, :);
    if ~isempty(CR)
        add = [cellstr(CR.File) cellstr(CR.SignalID) num2cell(CR.Photoperiod_h) ...
            cellstr(CR.Phase) cellstr(CR.BandName) cellstr(CR.Source) num2cell(CR.MeanBandPower_log10)];
        rows = [rows; add]; %#ok<AGROW>
    end
    UR = BCS(BCS.Source == SRC_UR & ismember(BCS.BandName, UR_BANDS), :);
    if ~isempty(UR)
        add = [cellstr(UR.File) cellstr(UR.SignalID) num2cell(UR.Photoperiod_h) ...
            cellstr(UR.Phase) cellstr(UR.BandName) cellstr(UR.Source) num2cell(UR.MeanBandPower_log10)];
        rows = [rows; add]; %#ok<AGROW>
    end
    TBP = cell2table(rows, 'VariableNames', hdr);
    TBP.File = string(TBP.File);
    TBP.SignalID = string(TBP.SignalID);
    TBP.Phase = string(TBP.Phase);
    TBP.BandName = string(TBP.BandName);
    TBP.Source = string(TBP.Source);
end

function result = run_lme_block_(BCS, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, cfg, LOG)
    result = struct('success', false, 'inferenceRaw', table(), 'inferenceFdr', table(), 'nSuccess', 0);
    if ~isfield(cfg, 'stats')
        return;
    end
    alpha = cfg.stats.alphaFdr;
    bandsAll = string(cfg.bands.allNames);

    inferenceRaw = table();
    nSuccess = 0;

    TpairL = Tpair;
    TpairL.Photoperiod_h = double(TpairL.Photoperiod_h);
    TpairL.SignalID = categorical(string(TpairL.SignalID));
    TpairL.File = categorical(string(TpairL.File));
    TpairL.Phase = categorical(string(TpairL.Phase));
    TpairL.UR_Band = categorical(string(TpairL.UR_Band));
    TpairL = add_sex_column_(TpairL);

    for b = 1:numel(UR_BANDS)
        ub = UR_BANDS(b);
        D = TpairL(TpairL.UR_Band == categorical(ub), :);
        D = drop_unknown_sex_(D);
        if height(D) < 3
            log_line_(LOG, 'Delta LME skipped %s: not enough rows', ub);
            continue;
        end
        fml = build_lme_formula_('Delta_log10', numel(categories(D.Phase)) > 1, has_usable_sex_(D));
        try
            mdl = fitlme(D, fml);
            inferenceRaw = append_inference_(inferenceRaw, mdl.Coefficients, ...
                "Delta_LME_Coefficients", "Delta", ub, "Delta_log10", fml, "Coefficients");
            inferenceRaw = append_inference_(inferenceRaw, anova(mdl), ...
                "Delta_LME_Anova", "Delta", ub, "Delta_log10", fml, "ANOVA");
            nSuccess = nSuccess + 1;
        catch ME
            log_line_(LOG, 'Delta LME failed %s: %s', ub, ME.message);
        end
    end

    TBP = build_bandpower_table_(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
    TBP.Photoperiod_h = double(TBP.Photoperiod_h);
    TBP.SignalID = categorical(string(TBP.SignalID));
    TBP.File = categorical(string(TBP.File));
    TBP.Phase = categorical(string(TBP.Phase));
    TBP.BandName = categorical(string(TBP.BandName));
    TBP = add_sex_column_(TBP);

    for b = 1:numel(bandsAll)
        bn = bandsAll(b);
        D = TBP(TBP.BandName == categorical(bn), :);
        D = drop_unknown_sex_(D);
        if height(D) < 3
            log_line_(LOG, 'Power LME skipped %s: not enough rows', bn);
            continue;
        end
        fml = build_lme_formula_('MeanBandPower_log10', numel(categories(D.Phase)) > 1, has_usable_sex_(D));
        try
            mdl = fitlme(D, fml);
            inferenceRaw = append_inference_(inferenceRaw, mdl.Coefficients, ...
                "Power_LME_Coefficients", "Power", bn, "MeanBandPower_log10", fml, "Coefficients");
            inferenceRaw = append_inference_(inferenceRaw, anova(mdl), ...
                "Power_LME_Anova", "Power", bn, "MeanBandPower_log10", fml, "ANOVA");
            nSuccess = nSuccess + 1;
        catch ME
            log_line_(LOG, 'Power LME failed %s: %s', bn, ME.message);
        end
    end

    result.nSuccess = nSuccess;
    result.inferenceRaw = inferenceRaw;
    if isempty(inferenceRaw)
        return;
    end

    inferenceFdr = apply_lme_bh_fdr_(inferenceRaw, alpha);
    result.inferenceFdr = inferenceFdr;
    result.success = nSuccess > 0;
end

function T = annotate_bcs_with_validation_(T, CarryForward, UR_BANDS, SRC_UR, applyAllPhase)
    T.Source = string(T.Source);
    T.BandName = string(T.BandName);
    T.File = string(T.File);
    T.SignalID = string(T.SignalID);
    T.Phase = string(T.Phase);
    if ismember('Photoperiod_h', T.Properties.VariableNames)
        T.Photoperiod_h = double(T.Photoperiod_h);
    end

    n = height(T);
    T.ValidatedRawUR = false(n, 1);
    T.ValidationStatus = repmat("NotApplicable", n, 1);
    isRawUR = T.Source == SRC_UR & ismember(T.BandName, UR_BANDS);
    T.ValidationStatus(isRawUR) = "RawUR_NotValidated";
    if isempty(CarryForward)
        return;
    end

    keysCarryAll = make_validation_keys_(CarryForward, false);
    for r = find(isRawUR).'
        keyAll = make_one_validation_key_(T.File(r), T.SignalID(r), T.Photoperiod_h(r), T.BandName(r), "");
        idx = find(keysCarryAll == keyAll, 1, 'first');
        if isempty(idx) || ~applyAllPhase
            keyExact = make_one_validation_key_(T.File(r), T.SignalID(r), T.Photoperiod_h(r), T.BandName(r), T.Phase(r));
            idx = find(make_validation_keys_(CarryForward, true) == keyExact, 1, 'first');
        end
        if ~isempty(idx)
            T.ValidatedRawUR(r) = true;
            T.ValidationStatus(r) = "CarryForward_ValidatedRawUR";
        end
    end
end

function Tused = filter_bcs_for_analysis_(T, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, requireValidatedUR)
    T.Source = string(T.Source);
    T.BandName = string(T.BandName);
    keepCR = T.Source == SRC_CR & T.BandName == CR_BAND;
    keepUR = T.Source == SRC_UR & ismember(T.BandName, UR_BANDS);
    if requireValidatedUR && ismember('ValidatedRawUR', T.Properties.VariableNames)
        keepUR = keepUR & T.ValidatedRawUR;
    end
    Tused = T(keepCR | keepUR, :);
end

function keys = make_validation_keys_(T, includePhase)
    n = height(T);
    keys = strings(n, 1);
    for i = 1:n
        phase = "";
        if includePhase
            phase = string(T.Phase(i));
        end
        pp = NaN;
        if ismember('Photoperiod_h', T.Properties.VariableNames)
            pp = double(T.Photoperiod_h(i));
        end
        keys(i) = make_one_validation_key_(string(T.File(i)), string(T.SignalID(i)), pp, string(T.BandName(i)), phase);
    end
end

function key = make_one_validation_key_(file, signalID, photoperiod_h, bandName, phase)
    key = lower(strjoin([string(file), string(signalID), sprintf('%.10g', double(photoperiod_h)), ...
        string(bandName), string(phase)], "|"));
end

function T = add_sex_column_(T)
    if ~ismember('SignalID', T.Properties.VariableNames)
        return;
    end
    sx = infer_sex_from_signalid_(string(T.SignalID));
    T.Sex = categorical(sx);
end

function D = drop_unknown_sex_(D)
    if ~ismember('Sex', D.Properties.VariableNames)
        return;
    end
    sx = string(D.Sex);
    D = D(~ismissing(sx) & sx ~= "Unknown", :);
end

function tf = has_usable_sex_(T)
    tf = false;
    if ~ismember('Sex', T.Properties.VariableNames)
        return;
    end
    sx = string(T.Sex);
    sx = sx(~ismissing(sx) & sx ~= "Unknown" & strlength(sx) > 0);
    tf = numel(unique(sx)) >= 2;
end

function fml = build_lme_formula_(responseName, hasPhase, hasSex)
    terms = "Photoperiod_h";
    if hasPhase
        terms = terms + " + Phase";
    end
    if hasSex
        terms = terms + " + Sex";
    end
    fml = char(string(responseName) + " ~ " + terms + " + (1|SignalID) + (1|File)");
end

function sx = infer_sex_from_signalid_(signalID)
    s = lower(string(signalID));
    sx = repmat("Unknown", numel(s), 1);
    isFemale = contains(s, "female") | contains(s, "_f_") | contains(s, "-f-") | ...
        startsWith(s, "f_") | startsWith(s, "f-") | endsWith(s, "_f") | endsWith(s, "-f");
    isMale = contains(s, "male") | contains(s, "_m_") | contains(s, "-m-") | ...
        startsWith(s, "m_") | startsWith(s, "m-") | endsWith(s, "_m") | endsWith(s, "-m");
    isMale = isMale & ~contains(s, "female");
    sx(isMale) = "Male";
    sx(isFemale) = "Female";
end

function Tinf = append_inference_(Tinf, T, fdrFamily, metricClass, bandName, responseName, formulaStr, tableType)
    T = coerce_lme_table_(T);
    if isempty(T) || height(T) == 0
        return;
    end
    n = height(T);
    term = strings(n, 1);
    if ismember('Name', T.Properties.VariableNames)
        term = string(T.Name);
    elseif ismember('Term', T.Properties.VariableNames)
        term = string(T.Term);
    else
        for i = 1:n
            term(i) = "Row_" + string(i);
        end
    end
    pRaw = nan(n, 1);
    for nm = {'pValue', 'PValue', 'Prob_F', 'ProbF'}
        if ismember(nm{1}, T.Properties.VariableNames)
            pRaw = double(T.(nm{1}));
            break;
        end
    end
    isIntercept = contains(lower(term), 'intercept');
    include = isfinite(pRaw) & ~isIntercept;
    if string(tableType) == "ANOVA"
        include = isfinite(pRaw);
    end

    block = table();
    block.FDRFamily = repmat(string(fdrFamily), n, 1);
    block.MetricClass = repmat(string(metricClass), n, 1);
    block.BandName = repmat(string(bandName), n, 1);
    block.Response = repmat(string(responseName), n, 1);
    block.TableType = repmat(string(tableType), n, 1);
    block.Term = term(:);
    block.Formula = repmat(string(formulaStr), n, 1);
    block.p_raw = pRaw(:);
    block.IncludeInFDR = include(:);

    if isempty(Tinf)
        Tinf = block;
    else
        Tinf = [Tinf; block]; %#ok<AGROW>
    end
end

function Tout = apply_lme_bh_fdr_(Tin, alphaVal)
    Tout = Tin;
    n = height(Tout);
    Tout.p_BH = nan(n, 1);
    Tout.Significant_BH = false(n, 1);
    fams = unique(string(Tout.FDRFamily));
    for f = 1:numel(fams)
        idxFam = string(Tout.FDRFamily) == fams(f) & Tout.IncludeInFDR & isfinite(Tout.p_raw);
        ii = find(idxFam);
        if isempty(ii), continue; end
        p = Tout.p_raw(ii);
        [ps, ord] = sort(p, 'ascend');
        m = numel(ps);
        ranks = (1:m)';
        qs = ps .* m ./ ranks;
        qs = flipud(cummin(flipud(qs)));
        qs(qs > 1) = 1;
        q = nan(m, 1);
        q(ord) = qs;
        Tout.p_BH(ii) = q;
        Tout.Significant_BH(ii) = q <= alphaVal;
    end
end

function T = coerce_lme_table_(T)
    if istable(T)
        return;
    end
    try
        T = table(T);
    catch
        T = table();
    end
end
