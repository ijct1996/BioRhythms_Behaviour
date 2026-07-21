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
    Tpair = add_sex_column_(Tpair);
    Tabs = build_absolute_summary_(BCS_used, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
    writetable(Tpair, xlsxOut, 'Sheet', 'CR_UR_Pairs_PerMouse');
    writetable(TpairSummary, xlsxOut, 'Sheet', 'CR_UR_Pairs_Summary');
    writetable(Tabs, xlsxOut, 'Sheet', 'AbsolutePower_Summary');

    sexResult = run_sex_descriptive_block_(BCS_used, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, LOG);
    writetable(sexResult.assignment, xlsxOut, 'Sheet', 'Sex_Assignment');
    writetable(sexResult.balance, xlsxOut, 'Sheet', 'Sex_Balance');
    writetable(sexResult.pairSummaryBySex, xlsxOut, 'Sheet', 'CR_UR_Pairs_Summary_BySex');
    writetable(sexResult.absSummaryBySex, xlsxOut, 'Sheet', 'AbsolutePower_Summary_BySex');
    fprintf('Sex-labelled mice: %d | unknown sex dropped from sex tables: %d\n', ...
        height(sexResult.assignment), sexResult.nUnknownSexRows);
    log_line_(LOG, 'Sex assignment rows: %d | usable for sex pass: %s', ...
        height(sexResult.assignment), string(sexResult.hasUsableSex));

    lmeResult = run_lme_block_(BCS_used, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, cfg, LOG, "AdditiveOnly");
    sexLmeResult = struct('success', false, 'inferenceRaw', table(), 'inferenceFdr', table(), 'nSuccess', 0);
    if sexResult.hasUsableSex && isfield(cfg, 'stats') && isfield(cfg.stats, 'sexInteractionLME') ...
            && cfg.stats.sexInteractionLME
        sexLmeResult = run_lme_block_(BCS_used, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, cfg, LOG, "PhotoperiodBySex");
        log_line_(LOG, 'Sex-interaction LME models: %d', sexLmeResult.nSuccess);
    end

    if lmeResult.success || sexLmeResult.success
        lmeXlsx = fullfile(dirs.Tables, 'AcrossPhotoperiod_LME_Outputs.xlsx');
        if isfile(lmeXlsx), delete(lmeXlsx); end
        inferenceRawAll = [lmeResult.inferenceRaw; sexLmeResult.inferenceRaw]; %#ok<AGROW>
        inferenceFdrAll = table();
        if ~isempty(lmeResult.inferenceFdr)
            inferenceFdrAll = lmeResult.inferenceFdr;
        end
        if ~isempty(sexLmeResult.inferenceFdr)
            inferenceFdrAll = [inferenceFdrAll; sexLmeResult.inferenceFdr]; %#ok<AGROW>
        end
        writetable(inferenceRawAll, lmeXlsx, 'Sheet', 'LME_Inference_Raw');
        writetable(inferenceFdrAll, lmeXlsx, 'Sheet', 'LME_Inference_BH_FDR');
        write_lme_fdr_subset_(lmeResult.inferenceFdr, lmeXlsx, "Delta", "ANOVA", 'LME_Anova_Delta_BH_FDR');
        write_lme_fdr_subset_(lmeResult.inferenceFdr, lmeXlsx, "Power", "ANOVA", 'LME_Anova_Power_BH_FDR');
        write_lme_fdr_subset_(lmeResult.inferenceFdr, lmeXlsx, "Delta", "Coefficients", 'LME_Coef_Delta_BH_FDR');
        write_lme_fdr_subset_(lmeResult.inferenceFdr, lmeXlsx, "Power", "Coefficients", 'LME_Coef_Power_BH_FDR');
        write_lme_fdr_family_subset_(sexLmeResult.inferenceFdr, lmeXlsx, "Delta_LME_SexInt", 'LME_SexInt_Delta_BH_FDR');
        write_lme_fdr_family_subset_(sexLmeResult.inferenceFdr, lmeXlsx, "Power_LME_SexInt", 'LME_SexInt_Power_BH_FDR');
        fprintf('LME workbook: %s\n', lmeXlsx);
        log_line_(LOG, 'LME workbook: %s', lmeXlsx);
        log_line_(LOG, 'LME inference rows: %d | BH-significant: %d', ...
            height(inferenceFdrAll), sum(inferenceFdrAll.Significant_BH));
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
    out.sex = sexResult;
    out.sexLme = sexLmeResult;
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

function result = run_lme_block_(BCS, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, cfg, LOG, sexMode)
    result = struct('success', false, 'inferenceRaw', table(), 'inferenceFdr', table(), 'nSuccess', 0);
    if nargin < 9 || isempty(sexMode)
        sexMode = "AdditiveOnly";
    end
    if ~isfield(cfg, 'stats')
        return;
    end
    alpha = cfg.stats.alphaFdr;
    bandsAll = string(cfg.bands.allNames);
    sexMode = string(sexMode);
    isSexInt = sexMode == "PhotoperiodBySex";
    coefFamilyPrefix = ternary_(isSexInt, "SexInt", "");
    anovaFamilyPrefix = coefFamilyPrefix;

    inferenceRaw = table();
    nSuccess = 0;

    TpairL = Tpair;
    TpairL.Photoperiod_h = double(TpairL.Photoperiod_h);
    TpairL.SignalID = categorical(string(TpairL.SignalID));
    TpairL.File = categorical(string(TpairL.File));
    TpairL.Phase = categorical(string(TpairL.Phase));
    TpairL.UR_Band = categorical(string(TpairL.UR_Band));
    if ~ismember('Sex', TpairL.Properties.VariableNames)
        TpairL = add_sex_column_(TpairL);
    end

    for b = 1:numel(UR_BANDS)
        ub = UR_BANDS(b);
        D = TpairL(TpairL.UR_Band == categorical(ub), :);
        D = drop_unknown_sex_(D);
        if height(D) < 3
            log_line_(LOG, '%s Delta LME skipped %s: not enough rows', sexMode, ub);
            continue;
        end
        hasSex = has_usable_sex_(D);
        if isSexInt && ~hasSex
            log_line_(LOG, 'Sex-interaction Delta LME skipped %s: need both sexes', ub);
            continue;
        end
        fml = build_lme_formula_('Delta_log10', numel(categories(D.Phase)) > 1, hasSex, sexMode);
        try
            mdl = fitlme(D, fml);
            inferenceRaw = append_inference_(inferenceRaw, mdl.Coefficients, ...
                build_lme_family_("Delta", "Coefficients", coefFamilyPrefix), "Delta", ub, "Delta_log10", fml, "Coefficients", sexMode);
            inferenceRaw = append_inference_(inferenceRaw, anova(mdl), ...
                build_lme_family_("Delta", "ANOVA", anovaFamilyPrefix), "Delta", ub, "Delta_log10", fml, "ANOVA", sexMode);
            nSuccess = nSuccess + 1;
        catch ME
            log_line_(LOG, '%s Delta LME failed %s: %s', sexMode, ub, ME.message);
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
            log_line_(LOG, '%s Power LME skipped %s: not enough rows', sexMode, bn);
            continue;
        end
        hasSex = has_usable_sex_(D);
        if isSexInt && ~hasSex
            log_line_(LOG, 'Sex-interaction Power LME skipped %s: need both sexes', bn);
            continue;
        end
        fml = build_lme_formula_('MeanBandPower_log10', numel(categories(D.Phase)) > 1, hasSex, sexMode);
        try
            mdl = fitlme(D, fml);
            inferenceRaw = append_inference_(inferenceRaw, mdl.Coefficients, ...
                build_lme_family_("Power", "Coefficients", coefFamilyPrefix), "Power", bn, "MeanBandPower_log10", fml, "Coefficients", sexMode);
            inferenceRaw = append_inference_(inferenceRaw, anova(mdl), ...
                build_lme_family_("Power", "ANOVA", anovaFamilyPrefix), "Power", bn, "MeanBandPower_log10", fml, "ANOVA", sexMode);
            nSuccess = nSuccess + 1;
        catch ME
            log_line_(LOG, '%s Power LME failed %s: %s', sexMode, bn, ME.message);
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

function result = run_sex_descriptive_block_(BCS, Tpair, CR_BAND, UR_BANDS, SRC_CR, SRC_UR, LOG)
    result = struct('assignment', table(), 'balance', table(), ...
        'pairSummaryBySex', table(), 'absSummaryBySex', table(), ...
        'hasUsableSex', false, 'nUnknownSexRows', 0);

    if ~ismember('Sex', Tpair.Properties.VariableNames)
        Tpair = add_sex_column_(Tpair);
    end
    sxAll = string(Tpair.Sex);
    result.nUnknownSexRows = sum(ismissing(sxAll) | sxAll == "Unknown");

    ids = unique(string(Tpair.SignalID));
    assignSex = strings(numel(ids), 1);
    for i = 1:numel(ids)
        rows = string(Tpair.SignalID) == ids(i);
        u = unique(sxAll(rows));
        u = u(~ismissing(u) & u ~= "Unknown");
        if numel(u) == 1
            assignSex(i) = u(1);
        else
            assignSex(i) = "Unknown";
        end
    end
    result.assignment = table(ids, assignSex, 'VariableNames', {'SignalID', 'Sex'});
    result.assignment = sortrows(result.assignment, {'Sex', 'SignalID'});

    known = result.assignment(~ismissing(result.assignment.Sex) & result.assignment.Sex ~= "Unknown", :);
    result.hasUsableSex = height(known) >= 2 && numel(unique(known.Sex)) >= 2;
    if ~result.hasUsableSex
        log_line_(LOG, 'Sex descriptive pass: insufficient labelled mice (nKnown=%d)', height(known));
        result.balance = table();
        result.pairSummaryBySex = table();
        result.absSummaryBySex = table();
        return;
    end

    Tknown = Tpair(~ismissing(sxAll) & sxAll ~= "Unknown", :);
    result.balance = build_sex_balance_table_(Tknown);
    result.pairSummaryBySex = build_cr_ur_pairs_summary_by_sex_(Tknown);
    result.absSummaryBySex = build_absolute_summary_by_sex_(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR);
    log_line_(LOG, 'Sex balance rows: %d | pair-by-sex rows: %d', ...
        height(result.balance), height(result.pairSummaryBySex));
end

function Tbal = build_sex_balance_table_(Tpair)
    T = Tpair(:, {'SignalID', 'Sex', 'Photoperiod_h', 'Phase', 'UR_Band'});
    T.SignalID = string(T.SignalID);
    T.Sex = string(T.Sex);
    T.Phase = string(T.Phase);
    T.UR_Band = string(T.UR_Band);
    G = findgroups(T.Sex, T.Photoperiod_h, T.Phase, T.UR_Band);
    sx = splitapply(@(x) x(1), T.Sex, G);
    pp = splitapply(@(x) x(1), T.Photoperiod_h, G);
    ph = splitapply(@(x) x(1), T.Phase, G);
    ub = splitapply(@(x) x(1), T.UR_Band, G);
    nMouse = splitapply(@(x) numel(unique(x)), T.SignalID, G);
    nObs = splitapply(@numel, T.SignalID, G);
    Tbal = table(sx, pp, ph, ub, nMouse, nObs, ...
        'VariableNames', {'Sex', 'Photoperiod_h', 'Phase', 'UR_Band', 'N_Mice', 'N_Obs'});
    Tbal = sortrows(Tbal, {'Sex', 'Phase', 'UR_Band', 'Photoperiod_h'});
end

function Tsum = build_cr_ur_pairs_summary_by_sex_(Tpair)
    G = findgroups(Tpair.Photoperiod_h, Tpair.Phase, Tpair.UR_Band, Tpair.Sex);
    meanCR = splitapply(@(x) mean(x, 'omitnan'), Tpair.CR_Log10, G);
    sdCR = splitapply(@(x) std(x, 'omitnan'), Tpair.CR_Log10, G);
    meanUR = splitapply(@(x) mean(x, 'omitnan'), Tpair.UR_Log10, G);
    sdUR = splitapply(@(x) std(x, 'omitnan'), Tpair.UR_Log10, G);
    meanD = splitapply(@(x) mean(x, 'omitnan'), Tpair.Delta_log10, G);
    sdD = splitapply(@(x) std(x, 'omitnan'), Tpair.Delta_log10, G);
    nMouse = splitapply(@(x) numel(unique(x)), string(Tpair.SignalID), G);
    [pp, ph, ub, sx] = splitapply(@(a, b, c, d) deal(a(1), b(1), c(1), d(1)), ...
        Tpair.Photoperiod_h, Tpair.Phase, Tpair.UR_Band, Tpair.Sex, G);
    Tsum = table(pp, ph, ub, sx, nMouse, meanCR, sdCR, meanUR, sdUR, meanD, sdD, ...
        'VariableNames', {'Photoperiod_h', 'Phase', 'UR_Band', 'Sex', 'N_Mice', ...
        'Mean_CR_Log10', 'SD_CR_Log10', 'Mean_UR_Log10', 'SD_UR_Log10', ...
        'Mean_Delta_log10', 'SD_Delta_log10'});
    Tsum.Sex = string(Tsum.Sex);
    Tsum = sortrows(Tsum, {'Sex', 'Phase', 'UR_Band', 'Photoperiod_h'});
end

function Tabs = build_absolute_summary_by_sex_(BCS, CR_BAND, UR_BANDS, SRC_CR, SRC_UR)
    BCS = add_sex_column_(BCS);
    sx = string(BCS.Sex);
    BCS = BCS(~ismissing(sx) & sx ~= "Unknown", :);
    bandsAll = [CR_BAND UR_BANDS];
    rows = cell(0, 6);
    for b = 1:numel(bandsAll)
        bn = bandsAll(b);
        if bn == CR_BAND
            D = BCS(BCS.Source == SRC_CR & BCS.BandName == bn, :);
        else
            D = BCS(BCS.Source == SRC_UR & BCS.BandName == bn, :);
        end
        if isempty(D), continue; end
        G = findgroups(D.Photoperiod_h, D.Phase, D.Sex);
        pp = splitapply(@(x) x(1), D.Photoperiod_h, G);
        ph = splitapply(@(x) x(1), D.Phase, G);
        sxG = splitapply(@(x) x(1), string(D.Sex), G);
        m = splitapply(@(x) mean(x, 'omitnan'), D.MeanBandPower_log10, G);
        s = splitapply(@(x) std(x, 'omitnan'), D.MeanBandPower_log10, G);
        n = numel(m);
        add = [num2cell(pp(:)) cellstr(ph(:)) cellstr(repmat(string(bn), n, 1)) ...
            cellstr(sxG(:)) num2cell(m(:)) num2cell(s(:))];
        rows = [rows; add]; %#ok<AGROW>
    end
    if isempty(rows)
        Tabs = table();
        return;
    end
    Tabs = cell2table(rows, 'VariableNames', ...
        {'Photoperiod_h', 'Phase', 'BandName', 'Sex', 'Mean_Log10', 'SD_Log10'});
    Tabs.Phase = string(Tabs.Phase);
    Tabs.BandName = string(Tabs.BandName);
    Tabs.Sex = string(Tabs.Sex);
    Tabs = sortrows(Tabs, {'Sex', 'Phase', 'BandName', 'Photoperiod_h'});
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

function fml = build_lme_formula_(responseName, hasPhase, hasSex, sexMode)
    if nargin < 4 || isempty(sexMode)
        sexMode = "AdditiveOnly";
    end
    terms = "Photoperiod_h";
    if hasPhase
        terms = terms + " + Phase";
    end
    if hasSex
        if string(sexMode) == "PhotoperiodBySex"
            terms = terms + " + Photoperiod_h:Sex + Sex";
        else
            terms = terms + " + Sex";
        end
    end
    fml = char(string(responseName) + " ~ " + terms + " + (1|SignalID) + (1|File)");
end

function fam = build_lme_family_(metricClass, tableType, prefix)
    if nargin < 3 || strlength(string(prefix)) == 0
        if string(tableType) == "ANOVA"
            fam = string(metricClass) + "_LME_Anova";
        else
            fam = string(metricClass) + "_LME_Coefficients";
        end
        return;
    end
    if string(tableType) == "ANOVA"
        fam = string(metricClass) + "_LME_" + string(prefix) + "_Anova";
    else
        fam = string(metricClass) + "_LME_" + string(prefix) + "_Coefficients";
    end
end

function y = ternary_(tf, a, b)
    if tf
        y = a;
    else
        y = b;
    end
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

function Tinf = append_inference_(Tinf, T, fdrFamily, metricClass, bandName, responseName, formulaStr, tableType, sexMode)
    if nargin < 9
        sexMode = "AdditiveOnly";
    end
    T = coerce_lme_table_(T);
    if isempty(T) || height(T) == 0
        return;
    end

    n = height(T);
    term = get_lme_terms_(T);
    pRaw = get_numeric_var_or_nan_(T, {'pValue', 'PValue', 'pVal', 'Prob_F', 'ProbF'});
    isIntercept = contains(lower(term), 'intercept');
    include = isfinite(pRaw) & ~isIntercept;
    if string(tableType) == "ANOVA"
        include = isfinite(pRaw);
    end
    if string(sexMode) == "PhotoperiodBySex" && string(tableType) == "Coefficients"
        isSexTerm = contains(lower(term), "sex") | contains(term, ":");
        include = include & isSexTerm;
    end
    if string(sexMode) == "PhotoperiodBySex" && string(tableType) == "ANOVA"
        isSexTerm = contains(lower(term), "sex") | contains(term, ":");
        include = include & isSexTerm;
    end

    block = table();
    block.FDRFamily = repmat(string(fdrFamily), n, 1);
    block.MetricClass = repmat(string(metricClass), n, 1);
    block.BandName = repmat(string(bandName), n, 1);
    block.Response = repmat(string(responseName), n, 1);
    block.TableType = repmat(string(tableType), n, 1);
    block.SexModel = repmat(string(sexMode), n, 1);
    block.Term = term(:);
    block.Formula = repmat(string(formulaStr), n, 1);
    block.p_raw = pRaw(:);
    block.IncludeInFDR = include(:);
    block.Estimate = get_numeric_var_or_nan_(T, {'Estimate'});
    block.SE = get_numeric_var_or_nan_(T, {'SE'});
    block.tStat = get_numeric_var_or_nan_(T, {'tStat'});
    block.FStat = get_numeric_var_or_nan_(T, {'FStat'});
    block.DF = get_numeric_var_or_nan_(T, {'DF'});
    block.DF1 = get_numeric_var_or_nan_(T, {'DF1'});
    block.DF2 = get_numeric_var_or_nan_(T, {'DF2'});

    if isempty(Tinf)
        Tinf = block;
    else
        Tinf = [Tinf; block]; %#ok<AGROW>
    end
end

function term = get_lme_terms_(T)
    n = height(T);
    term = strings(n, 1);
    if ismember('Name', T.Properties.VariableNames)
        term = string(T.Name);
        return;
    end
    if ismember('Term', T.Properties.VariableNames)
        term = string(T.Term);
        return;
    end
    if ~isempty(T.Properties.RowNames) && numel(T.Properties.RowNames) == n
        term = string(T.Properties.RowNames(:));
        return;
    end
    for i = 1:n
        term(i) = "Row_" + string(i);
    end
end

function x = get_numeric_var_or_nan_(T, nameCandidates)
    n = height(T);
    x = nan(n, 1);
    for k = 1:numel(nameCandidates)
        nm = char(nameCandidates{k});
        if ismember(nm, T.Properties.VariableNames)
            v = T.(nm);
            if iscell(v)
                v = string(v);
            end
            x = double(v);
            x = x(:);
            if numel(x) ~= n
                x = nan(n, 1);
            end
            return;
        end
    end
end

function write_lme_fdr_family_subset_(inferenceFdr, xlsxPath, familyPrefix, sheetName)
    if isempty(inferenceFdr)
        writecell({'No rows.'}, xlsxPath, 'Sheet', sheetName);
        return;
    end
    fams = string(inferenceFdr.FDRFamily);
    subset = inferenceFdr(startsWith(fams, string(familyPrefix)), :);
    if isempty(subset)
        writecell({'No rows.'}, xlsxPath, 'Sheet', sheetName);
    else
        writetable(subset, xlsxPath, 'Sheet', sheetName);
    end
end

function write_lme_fdr_subset_(inferenceFdr, xlsxPath, metricClass, tableType, sheetName)
    if isempty(inferenceFdr)
        return;
    end
  subset = inferenceFdr(string(inferenceFdr.MetricClass) == string(metricClass) & ...
        string(inferenceFdr.TableType) == string(tableType), :);
    if isempty(subset)
        writecell({'No rows.'}, xlsxPath, 'Sheet', sheetName);
    else
        writetable(subset, xlsxPath, 'Sheet', sheetName);
    end
end

function Tout = apply_lme_bh_fdr_(Tin, alphaVal)
    Tout = Tin;
    n = height(Tout);
    Tout.p_BH = nan(n, 1);
    Tout.Significant_BH = false(n, 1);
    Tout.FDR_Rank = nan(n, 1);
    Tout.FDR_m = nan(n, 1);
    Tout.FDR_Alpha = repmat(alphaVal, n, 1);
    fams = unique(string(Tout.FDRFamily));
    fams = fams(~ismissing(fams));
    for f = 1:numel(fams)
        idxFam = string(Tout.FDRFamily) == fams(f) & Tout.IncludeInFDR & isfinite(Tout.p_raw);
        ii = find(idxFam);
        if isempty(ii), continue; end
        p = Tout.p_raw(ii);
        [ps, ord] = sort(p, 'ascend');
        m = numel(ps);
        ranks = (1:m)';
        qs = ps .* m ./ ranks;
        for i = m-1:-1:1
            qs(i) = min(qs(i), qs(i + 1));
        end
        qs(qs > 1) = 1;
        q = nan(m, 1);
        q(ord) = qs;
        rankOut = nan(m, 1);
        rankOut(ord) = ranks;
        Tout.p_BH(ii) = q;
        Tout.Significant_BH(ii) = q <= alphaVal;
        Tout.FDR_Rank(ii) = rankOut;
        Tout.FDR_m(ii) = m;
    end
end

function T = coerce_lme_table_(T)
    if istable(T)
        return;
    end
    % fitlme / anova return titleddataset in R2025b; table(T) collapses to a
    % single unusable column. dataset2table preserves Name/Term and pValue.
    try
        T = dataset2table(T);
        return;
    catch
    end
    try
        if isa(T, 'dataset')
            T = dataset2table(T);
            return;
        end
    catch
    end
    T = table();
end
