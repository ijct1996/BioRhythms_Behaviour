function indexTable = extended_script10_export_supplemental_tables(paths, destDirs, cfg, LOG)
%EXTENDED_SCRIPT10_EXPORT_SUPPLEMENTAL_TABLES Copy/export all manuscript supp tables.
%
%   Writes full workbooks + key CSV extracts to:
%     destDirs.packageTables   (Package_Manuscript/Tables)
%     destDirs.script10Tables  (Script10_*/Tables/Supplemental) — when set
%
%   indexTable = extended_script10_export_supplemental_tables(paths, destDirs, cfg, LOG)
%
%   Read-only aggregation from Scripts 4–7 (+ optional Script 9). Does not recompute stats.
%   Script 11 outputs stay in Script11_DominantPeriod_* (not bundled here).

    if nargin < 4, LOG = -1; end
    if nargin < 3 || isempty(cfg)
        cfg = extended_defaults();
    end
    exportCsv = true;
    if isfield(cfg, 'script10') && isfield(cfg.script10, 'exportSupplementalTables')
        if ~cfg.script10.exportSupplementalTables
            indexTable = table();
            return;
        end
    end
    if isfield(cfg, 'script10') && isfield(cfg.script10, 'exportSupplementalCsv')
        exportCsv = logical(cfg.script10.exportSupplementalCsv);
    end

    pkgTab = destDirs.packageTables;
    wbDir = fullfile(pkgTab, 'Workbooks');
    csvDir = fullfile(pkgTab, 'CSV');
    extended_period_gate_ensure_dir(wbDir);
    extended_period_gate_ensure_dir(csvDir);

    script10Tab = '';
    if isfield(destDirs, 'script10Tables') && strlength(string(destDirs.script10Tables)) > 0
        script10Tab = char(destDirs.script10Tables);
        extended_period_gate_ensure_dir(script10Tab);
        extended_period_gate_ensure_dir(fullfile(script10Tab, 'Workbooks'));
        extended_period_gate_ensure_dir(fullfile(script10Tab, 'CSV'));
    end

    optPaths = script10_supp_resolve_optional_paths_(paths.cohortRoot, cfg);
    catalog = script10_supp_table_catalog_(paths, optPaths);
    indexRows = {};

    script10_supp_log_(LOG, 'Exporting %d supplemental table entries', numel(catalog));

    for i = 1:numel(catalog)
        entry = catalog(i);
        srcXlsx = char(entry.SourceXlsx);
        if ~isfile(srcXlsx)
            indexRows(end + 1, :) = {string(entry.Category), string(entry.Sheet), ...
                string(entry.CsvName), "missing_source", string(entry.ManuscriptUse), ...
                string(srcXlsx), ""}; %#ok<AGROW>
            script10_supp_log_(LOG, 'SKIP (missing): %s / %s', basename_(srcXlsx), char(entry.Sheet));
            continue;
        end

        % Workbook copy (once per source file)
        wbName = char(entry.WorkbookCopyName);
        wbDestPkg = fullfile(wbDir, wbName);
        if ~isfile(wbDestPkg)
            script10_supp_copy_file_(srcXlsx, wbDestPkg, LOG);
        end
        if strlength(string(script10Tab)) > 0
            wbDestS10 = fullfile(script10Tab, 'Workbooks', wbName);
            if ~isfile(wbDestS10)
                script10_supp_copy_file_(srcXlsx, wbDestS10, LOG);
            end
        end

        if ~entry.ExportCsv
            indexRows(end + 1, :) = {string(entry.Category), string(entry.Sheet), ...
                "", "workbook_only", string(entry.ManuscriptUse), string(srcXlsx), wbName}; %#ok<AGROW>
            continue;
        end

        csvName = char(entry.CsvName);
        if strlength(string(csvName)) == 0
            csvName = sprintf('%s__%s.csv', erase(wbName, '.xlsx'), char(entry.Sheet));
        end
        csvPathPkg = fullfile(csvDir, csvName);
        status = script10_supp_export_sheet_csv_(srcXlsx, char(entry.Sheet), csvPathPkg, LOG);

        csvPathS10 = '';
        if strlength(string(script10Tab)) > 0 && isfile(csvPathPkg)
            csvPathS10 = fullfile(script10Tab, 'CSV', csvName);
            script10_supp_copy_file_(csvPathPkg, csvPathS10, LOG);
        end

        indexRows(end + 1, :) = {string(entry.Category), string(entry.Sheet), ...
            string(csvName), string(status), string(entry.ManuscriptUse), ...
            string(srcXlsx), wbName}; %#ok<AGROW>
    end

    indexHdr = {'Category','Sheet','CsvFile','Status','ManuscriptUse','SourceXlsx','WorkbookCopy'};
    indexTable = script10_supp_rows_to_table_(indexRows, indexHdr);
    if ~isempty(indexTable)
        writetable(indexTable, fullfile(pkgTab, 'SUPPLEMENTAL_TABLE_INDEX.csv'));
        if strlength(string(script10Tab)) > 0
            writetable(indexTable, fullfile(script10Tab, 'SUPPLEMENTAL_TABLE_INDEX.csv'));
        end
    end

    script10_supp_write_readme_(pkgTab, indexTable);
    if strlength(string(script10Tab)) > 0
        script10_supp_write_readme_(script10Tab, indexTable);
    end

    script10_supp_log_(LOG, 'Supplemental tables exported to %s (%d index rows)', pkgTab, height(indexTable));
end

%% Catalog
function catalog = script10_supp_table_catalog_(paths, optPaths)
    catalog = struct('Category', {}, 'SourceXlsx', {}, 'Sheet', {}, ...
        'CsvName', {}, 'WorkbookCopyName', {}, 'ExportCsv', {}, 'ManuscriptUse', {});

    gate = paths.gateXlsx;
    lmeDesc = paths.lmeDescriptiveXlsx;
    lmeInf = paths.lmeInferenceXlsx;
    resync = paths.resyncXlsx;
    gradient = paths.resyncGradientXlsx;
    profiles = paths.profilesXlsx;
    profilesBase = basename_(profiles);

    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'CarryForward_Periods', ...
        'Script04_CarryForward_Periods.csv', 'HSubSupported_PeriodMap.xlsx', true, ...
        'Validated Raw UR candidates included downstream (±15% SEL_P360 CarryForward)')];
    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'Retention_ByBand', ...
        'Script04_Retention_ByBand.csv', 'HSubSupported_PeriodMap.xlsx', true, ...
        'CarryForward retention by band (Fig 2E citation / methods)')];
    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'Retention_ByPhotoperiod', ...
        'Script04_Retention_ByPhotoperiod.csv', 'HSubSupported_PeriodMap.xlsx', true, ...
        'CarryForward retention by photoperiod')];
    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'Retention_ByPhotoBand', ...
        'Script04_Retention_ByPhotoBand.csv', 'HSubSupported_PeriodMap.xlsx', true, ...
        'CarryForward retention by photoperiod × band')];
    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'Matched_Periods_All', ...
        'Script04_Matched_Periods_All.csv', 'HSubSupported_PeriodMap.xlsx', true, ...
        'Raw vs HSub period match audit (full ladder view)')];
    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'QC_Flags', ...
        'Script04_QC_Flags.csv', 'HSubSupported_PeriodMap.xlsx', true, ...
        'Harmonic-sensitive / QC-flagged CarryForward rows')];
    catalog = [catalog; script10_supp_entry_('Methods_QC', gate, 'AllPeriodCandidates', ...
        '', 'HSubSupported_PeriodMap.xlsx', false, ...
        'All ridge period candidates (workbook only; large)')];

    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'CR_UR_Pairs_Summary', ...
        'CR_UR_Pairs_Summary.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Fig 3A-B descriptive CR–UR co-expression means')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'CR_UR_Pairs_PerMouse', ...
        'CR_UR_Pairs_PerMouse.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Per-mouse CR–UR pair descriptives (supp audit)')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'AbsolutePower_Summary', ...
        'AbsolutePower_Summary.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Fig 3A absolute band-power summary')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'Sex_Balance', ...
        'Sex_Balance.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Sample size / sex balance (methods + Script 9 supp context)')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'Sex_Assignment', ...
        'Sex_Assignment.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Mouse sex assignment audit')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'CR_UR_Pairs_Summary_BySex', ...
        'CR_UR_Pairs_Summary_BySex.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Sex-stratified co-expression descriptives (Script 9 supp)')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'AbsolutePower_Summary_BySex', ...
        'AbsolutePower_Summary_BySex.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Sex-stratified absolute power (Script 9 supp)')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'BandCondSummary_AnalysisUsed', ...
        'BandConditionSummary_AnalysisUsed.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'Band×condition rows fed to Script 6 LME')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeDesc, 'ValidationMap_Used', ...
        'ValidationMap_Used.csv', 'AcrossPhotoperiod_Outputs.xlsx', true, ...
        'CarryForward rows used in Script 6')];

    catalog = [catalog; script10_supp_entry_('CoExpression', lmeInf, 'LME_Coef_Delta_BH_FDR', ...
        'LME_Coef_Delta_BH_FDR.csv', 'AcrossPhotoperiod_LME_Outputs.xlsx', true, ...
        'Fig 3C confirmatory LME — UR−CR balance (Δ) coefficients + BH-FDR')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeInf, 'LME_Coef_Power_BH_FDR', ...
        'LME_Coef_Power_BH_FDR.csv', 'AcrossPhotoperiod_LME_Outputs.xlsx', true, ...
        'Confirmatory LME — absolute UR power coefficients + BH-FDR')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeInf, 'LME_Inference_BH_FDR', ...
        'LME_Inference_BH_FDR.csv', 'AcrossPhotoperiod_LME_Outputs.xlsx', true, ...
        'Full Script 6 LME inference table (all terms, BH-FDR)')];
    catalog = [catalog; script10_supp_entry_('CoExpression', lmeInf, 'LME_Inference_Raw', ...
        'LME_Inference_Raw.csv', 'AcrossPhotoperiod_LME_Outputs.xlsx', true, ...
        'Script 6 LME raw p-values (pre-FDR audit)')];

    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'Resync_PrimaryStats_BH_FDR', ...
        'Resync_PrimaryStats_BH_FDR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Confirmatory transition resync (ΔR>0); stars on Supp LD Pre/Post R bars')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'Resync_PrimaryStats', ...
        'Resync_PrimaryStats.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Transition resync primary stats (pre-FDR)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'RealVsPseudoStats_BH_FDR', ...
        'Resync_RealVsPseudoStats_BH_FDR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Real vs pseudo-transition control (BH-FDR)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'DeltaR_Summary', ...
        'DeltaR_Summary.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Group-level pre/post phase-coherence R summaries')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'PrePostCoherence', ...
        'PrePostCoherence.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Transition-aligned pre/post coherence summaries')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'CandidateDeltaR', ...
        'CandidateDeltaR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Candidate-level ΔR (supp audit / LD bar plot source)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'CandidatePrePost', ...
        'CandidatePrePost.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Candidate-level pre/post coherence')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'BinnedCoherence', ...
        'BinnedCoherence.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Binned transition phase-coherence (methods / exploratory)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'RidgePowerStats_BH_FDR', ...
        'RidgePowerStats_BH_FDR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Secondary ridge-power pre/post family (BH-FDR)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'RidgePeriodStats_BH_FDR', ...
        'RidgePeriodStats_BH_FDR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Secondary ridge-period pre/post family (BH-FDR)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'ValidatedCandidates_Used', ...
        'Resync_ValidatedCandidates_Used.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'CarryForward candidate lookup used in Script 5')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'TransitionPhase_Long', ...
        '', 'Ultradian_RidgePhase_Resync_Output.xlsx', false, ...
        'Long transition-phase series (workbook only; large)')];

    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'LLProj_PrimaryStats_BH_FDR', ...
        'LLProj_PrimaryStats_BH_FDR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'LL projected-dark primary ΔR family (BH-FDR)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'L22_vs_LLProj_BH_FDR', ...
        'L22_vs_LLProjectedStats_BH_FDR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'Real L22 vs projected LL aftereffect (BH-FDR)')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'LLProj_DeltaR', ...
        'LLProj_DeltaR.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'LL projected group ΔR summaries')];
    catalog = [catalog; script10_supp_entry_('TransitionResync', resync, 'LL_NoTransitionSummary', ...
        'LL_NoTransitionSummary.csv', 'Ultradian_RidgePhase_Resync_Output.xlsx', true, ...
        'LL candidates excluded from true transition inference')];

    catalog = [catalog; script10_supp_entry_('TransitionResync', gradient, 'Summary', ...
        'TransitionEffect_vs_Photoperiod_Summary.csv', 'TransitionEffect_vs_Photoperiod.xlsx', true, ...
        'Photoperiod gradient of transition effects (exploratory summary)')];

    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'ClusterSummary', ...
        'ClusterSummary.csv', profilesBase, true, ...
        'Validated period clusters (main + supp cluster lock)')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'ClusterMembership', ...
        'ClusterMembership.csv', profilesBase, true, ...
        'Candidate → cluster map (n annotations on Figs 4–7)')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'AllCandidatePeriods', ...
        'AllCandidatePeriods.csv', profilesBase, true, ...
        'CarryForward candidate median ridge periods')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'PhaseCoherence_24h', ...
        'PhaseCoherence_24h.csv', profilesBase, true, ...
        '24h ZT phase-coherence bins (Figs 5/7 source data)')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'PhaseCoherence_DL_LD', ...
        'PhaseCoherence_DL_LD.csv', profilesBase, true, ...
        'Transition-aligned phase coherence (exploratory / legacy panels)')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'ActivityComponent_24h', ...
        'ActivityComponent_24h.csv', profilesBase, true, ...
        '24h ZT z-scored activity components (Figs 4/6 source data)')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'RidgePower_24h', ...
        'RidgePower_24h.csv', profilesBase, true, ...
        '24h ZT ridge-power profiles (supp / methods audit)')];
    catalog = [catalog; script10_supp_entry_('Profiles', profiles, 'ActivityWarnings', ...
        'ActivityWarnings.csv', profilesBase, true, ...
        'Activity filter warnings (QC audit)')];

    if isfield(optPaths, 'script9Manifest') && isfile(optPaths.script9Manifest)
        catalog = [catalog; script10_supp_entry_('SupplementaryFigures', optPaths.script9Manifest, 'Manifest', ...
            'Script09_Supplementary_FigureManifest.csv', basename_(optPaths.script9Manifest), true, ...
            'Script 9 sex/cluster supp figure index (when Script 9 has been run)')];
    end
end

function e = script10_supp_entry_(category, xlsxPath, sheet, csvName, wbCopyName, exportCsv, manuscriptUse)
    e = struct();
    e.Category = string(category);
    e.SourceXlsx = string(xlsxPath);
    e.Sheet = string(sheet);
    e.CsvName = string(csvName);
    e.WorkbookCopyName = string(wbCopyName);
    e.ExportCsv = logical(exportCsv);
    e.ManuscriptUse = string(manuscriptUse);
end

%% Optional paths (Script 9 / 11)
function opt = script10_supp_resolve_optional_paths_(cohortRoot, ~)
    opt = struct();
    opt.script9Root = script10_supp_find_dir_(cohortRoot, 'Script9_SupplementaryFigures_*');
    opt.script9Manifest = '';
    if strlength(string(opt.script9Root)) > 0
        hits = dir(fullfile(char(opt.script9Root), '**', '*Manifest*.xlsx'));
        if ~isempty(hits)
            opt.script9Manifest = fullfile(hits(1).folder, hits(1).name);
        end
    end
end

%% I/O helpers
function status = script10_supp_export_sheet_csv_(xlsxPath, sheetName, csvPath, LOG)
    status = "empty";
    try
        opts = detectImportOptions(xlsxPath, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
        T = readtable(xlsxPath, opts);
    catch ME
        status = "read_error";
        script10_supp_log_(LOG, 'WARN sheet %s/%s: %s', basename_(xlsxPath), sheetName, ME.message);
        return;
    end
    if isempty(T) || height(T) == 0
        return;
    end
    try
        writetable(T, csvPath);
        status = "ok";
        script10_supp_log_(LOG, 'CSV %s (%d rows)', basename_(csvPath), height(T));
    catch ME
        status = "write_error";
        script10_supp_log_(LOG, 'WARN write %s: %s', csvPath, ME.message);
    end
end

function script10_supp_copy_file_(src, dest, LOG)
    if ~isfile(src), return; end
    extended_period_gate_ensure_dir(fileparts(dest));
    copyfile(src, dest, 'f');
    script10_supp_log_(LOG, 'Copy %s -> %s', basename_(src), dest);
end

function script10_supp_write_readme_(tabDir, indexTable)
    fid = fopen(fullfile(tabDir, 'README_SUPPLEMENTAL_TABLES.md'), 'w', 'n', 'UTF-8');
    if fid < 0, return; end
    fprintf(fid, '# Supplemental tables — manuscript handoff\n\n');
    fprintf(fid, 'Auto-exported by **Extended Script 10** (`extended_script10_export_supplemental_tables`).\n\n');
    fprintf(fid, '## Layout\n\n');
    fprintf(fid, '- `Workbooks/` — full source `.xlsx` files from Scripts 4–7 (and Script 9/11 when present)\n');
    fprintf(fid, '- `CSV/` — lightweight extracts for submission / sharing\n');
    fprintf(fid, '- `SUPPLEMENTAL_TABLE_INDEX.csv` — row-level export log\n\n');
    fprintf(fid, '## Categories\n\n');
    fprintf(fid, '| Category | Source | Typical manuscript use |\n');
    fprintf(fid, '|----------|--------|------------------------|\n');
    fprintf(fid, '| **Methods_QC** | Script 4 CarryForward gate | Validation, retention, harmonic QC |\n');
    fprintf(fid, '| **CoExpression** | Script 6 LME | Fig 3 descriptives + confirmatory LME supp |\n');
    fprintf(fid, '| **TransitionResync** | Script 5 | Supp LD Pre/Post R stats, LL projected controls |\n');
    fprintf(fid, '| **Profiles** | Script 7 | Figs 4–7 source tables, cluster definitions |\n');
    fprintf(fid, '| **SupplementaryFigures** | Script 9 | Sex-stratified supp figure manifest |\n\n');
    fprintf(fid, '## Main-text vs supplement\n\n');
    fprintf(fid, '- **Main text cite (Script 10 package):** `LME_Coef_*_BH_FDR`, `Resync_PrimaryStats_BH_FDR`, `ClusterSummary`\n');
    fprintf(fid, '- **Dominant UR periods (Script 11, separate folder):** `Script11_DominantPeriod_*/Tables/DominantPeriod_*`\n');
    fprintf(fid, '- **Supplement tables:** remaining CSVs (full audit, sex-stratified descriptives, secondary FDR families)\n');
    fprintf(fid, '- **Workbook-only rows:** large long-format sheets (see index `Status = workbook_only`)\n\n');
    if ~isempty(indexTable) && height(indexTable) > 0
        fprintf(fid, '## Export summary\n\n');
        fprintf(fid, 'Total indexed exports: **%d**\n', height(indexTable));
    end
    fclose(fid);
end

function d = script10_supp_find_dir_(root, pattern)
    d = '';
    hits = dir(fullfile(root, pattern));
    hits = hits([hits.isdir]);
    if ~isempty(hits)
        d = fullfile(hits(1).folder, hits(1).name);
    end
end

function T = script10_supp_rows_to_table_(rows, hdr)
    if isempty(rows)
        T = cell2table(cell(0, numel(hdr)), 'VariableNames', hdr);
        return;
    end
    T = cell2table(rows, 'VariableNames', hdr);
end

function script10_supp_log_(LOG, fmt, varargin)
    line = sprintf(fmt, varargin{:});
    fprintf('%s\n', line);
    if ~isempty(LOG) && LOG > 0
        fprintf(LOG, '%s\n', line);
    end
end

function b = basename_(p)
    [~, b, ext] = fileparts(char(p));
    if ~isempty(ext)
        b = [b ext];
    end
end
