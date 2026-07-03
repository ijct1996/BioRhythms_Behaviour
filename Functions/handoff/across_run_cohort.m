function summary = across_run_cohort(handoffDir, outputFolder, cohortTag, plotMode)
%ACROSS_RUN_COHORT Script 3 — RAW summaries + HSub-validated UR (Handoff/ only).
    if nargin < 4, plotMode = 'development'; end
    theme = plot_config(plotMode);
    theme = plot_theme_ensure_scalogram(theme);
    cfg = core_defaults();

    cohortHandoff = handoffDir;
    if ~isfolder(cohortHandoff)
        error('across_run_cohort:NoHandoff', 'Handoff folder not found: %s', cohortHandoff);
    end

    outFolder = fullfile(outputFolder, ['03_AcrossPhotoperiod_' cohortTag]);
    figFolder = fullfile(outFolder, 'Figures');
    tblFolder = fullfile(outFolder, 'Tables');
    ensure_dir(figFolder);
    ensure_dir(tblFolder);

    entries = across_load_handoff_cohort(cohortHandoff);

    allPeriod = vertcat(entries.periodTable);
    allBand = vertcat(entries.bandPower);
    ratioTable = across_compute_coexpression(allBand, cfg);

    % HSub-validated ultradian layer (confirmatory): UR from residual, CR from RAW
    [hsubBand, hsubPeriod, hsubValLog] = across_compute_hsub_validated_ur(entries, cfg);
    hsubRatio = across_compute_coexpression(hsubBand, cfg);

    writetable(allPeriod, fullfile(tblFolder, 'PeriodPeaks_all_photoperiods_RAW.xlsx'));
    writetable(ratioTable, fullfile(tblFolder, 'Coexpression_ratios_RAW.xlsx'));
    writetable(hsubBand, fullfile(tblFolder, 'BandPower_hsub_validated.xlsx'));
    writetable(hsubPeriod, fullfile(tblFolder, 'PeriodPeaks_hsub_validated.xlsx'));
    writetable(hsubValLog, fullfile(tblFolder, 'PeriodPeaks_hsub_validation_log.xlsx'));
    writetable(hsubRatio, fullfile(tblFolder, 'Coexpression_ratios_hsub_validated.xlsx'));

    ext = theme.exportFormat;
    across_plot_coexpression(ratioTable, fullfile(figFolder, ['Coexpression_CR_UR_RAW.' ext]), theme, ...
        'titleStr', 'Co-expression RAW (exploratory): CR vs UR band power', ...
        'ylabelStr', 'CR_{20-28} / UR_{total}');
    across_plot_period_comparison(allPeriod, fullfile(figFolder, ['Period_comparison_RAW.' ext]), theme, ...
        'titleStr', 'Period comparison RAW (exploratory): rank-1 peak', ...
        'ylabelStr', 'Dominant period (h) — rank 1', 'ylimRange', [0 26]);

    across_plot_coexpression(hsubRatio, fullfile(figFolder, ['Coexpression_CR_UR_hsub_validated.' ext]), theme, ...
        'titleStr', 'Co-expression HSub-validated: RAW CR vs residual UR', ...
        'ylabelStr', 'CR_{20-28} (RAW) / UR_{total} (residual)');
    across_plot_period_comparison(hsubPeriod, fullfile(figFolder, ['Period_comparison_hsub_validated.' ext]), theme, ...
        'titleStr', 'HSub-validated ultradian periods (confirmatory): rank-1 UR peak', ...
        'ylabelStr', 'Validated UR period (h) — rank 1', 'ylimRange', [0 18]);

    distRawPaths = across_plot_period_distributions(allPeriod, figFolder, ...
        'Period_distribution_RAW', theme, ...
        'titlePrefix', 'RAW rank-1 period distribution', ...
        'ylabelStr', 'Dominant period (h) — rank 1', 'ylimRange', [0 26]);
    distHsubPaths = across_plot_period_distributions(hsubPeriod, figFolder, ...
        'Period_distribution_hsub_validated', theme, ...
        'titlePrefix', 'HSub-validated UR period distribution', ...
        'ylabelStr', 'Validated UR period (h) — rank 1', 'ylimRange', [0 18]);

    stitchedPaths = across_plot_stitched_scalograms(entries, figFolder, cohortTag, theme);
    stitchedHsubPaths = across_plot_stitched_hsub_scalograms(entries, figFolder, cohortTag, theme);

    summary = struct('cohort', cohortTag, 'nPhotoperiods', numel(entries), ...
        'outFolder', outFolder, 'ratioTableRAW', ratioTable, ...
        'ratioTableHsubValidated', hsubRatio, ...
        'periodTableHsubValidated', hsubPeriod, ...
        'periodDistributionRAW', {distRawPaths}, ...
        'periodDistributionHsub', {distHsubPaths}, ...
        'stitchedScalograms', {stitchedPaths}, ...
        'stitchedHsubScalograms', {stitchedHsubPaths});
    fprintf('Across-photoperiod complete: %s (%d photoperiods)\n', outFolder, numel(entries));
    fprintf('  HSub-validated UR rows: band=%d, period peaks=%d\n', ...
        height(hsubBand), height(hsubPeriod));
end
