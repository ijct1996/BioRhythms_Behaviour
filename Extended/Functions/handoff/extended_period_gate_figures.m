function extended_period_gate_figures(M, PrimaryHSub, figDir, dpi, ext, primaryMode)
%EXTENDED_PERIOD_GATE_FIGURES Validation figures (JPEG) for CarryForward gate.

    extended_period_gate_ensure_dir(figDir);

    % Figure 1: paired Raw vs primary HSub period scatter.
    hasPair = isfinite(M.RawPeriod_h) & isfinite(M.PrimaryHSubPeriod_h);
    if any(hasPair)
        f = figure('Color', 'w', 'Position', [100 100 850 700]);
        scatter(M.RawPeriod_h(hasPair), M.PrimaryHSubPeriod_h(hasPair), 36, 'filled');
        hold on;
        minP = min([M.RawPeriod_h(hasPair); M.PrimaryHSubPeriod_h(hasPair)], [], 'omitnan');
        maxP = max([M.RawPeriod_h(hasPair); M.PrimaryHSubPeriod_h(hasPair)], [], 'omitnan');
        x = linspace(max(0.1, minP * 0.9), maxP * 1.1, 100);
        plot(x, x, 'k-', 'LineWidth', 1.2);
        plot(x, x * 1.15, 'k--', 'LineWidth', 0.8);
        plot(x, x / 1.15, 'k--', 'LineWidth', 0.8);
        xlabel('Raw ridge period (h)', 'FontWeight', 'bold');
        ylabel(sprintf('%s ridge period (h)', char(primaryMode)), 'FontWeight', 'bold');
        title('Raw versus Selective-HSub ultradian period agreement', 'FontWeight', 'bold');
        set(gca, 'FontName', 'Times New Roman', 'Box', 'off', 'TickDir', 'out');
        axis square;
        exportgraphics(f, fullfile(figDir, ['Raw_vs_' char(primaryMode) '_PeriodScatter' ext]), ...
            'Resolution', dpi);
        close(f);
    end

    % Figure 2: histogram overlay.
    rawP = M.RawPeriod_h(M.RawPassQC & isfinite(M.RawPeriod_h));
    hsubP = PrimaryHSub.MedianRidgePeriod_h( ...
        extended_period_gate_passes_qc(PrimaryHSub, 0, 0, false) & ...
        isfinite(PrimaryHSub.MedianRidgePeriod_h));
    if ~isempty(rawP) || ~isempty(hsubP)
        f = figure('Color', 'w', 'Position', [100 100 900 650]);
        hold on;
        edges = 1:0.5:18;
        if ~isempty(rawP)
            histogram(rawP, edges, 'Normalization', 'probability', ...
                'DisplayStyle', 'stairs', 'LineWidth', 1.5);
        end
        if ~isempty(hsubP)
            histogram(hsubP, edges, 'Normalization', 'probability', ...
                'DisplayStyle', 'stairs', 'LineWidth', 1.5);
        end
        xlabel('Period (h)', 'FontWeight', 'bold');
        ylabel('Proportion of candidates', 'FontWeight', 'bold');
        title('Raw and Selective-HSub candidate period distributions', 'FontWeight', 'bold');
        legend({'Raw', 'Selective-HSub'}, 'Location', 'best', 'Box', 'off');
        set(gca, 'FontName', 'Times New Roman', 'Box', 'off', 'TickDir', 'out');
        exportgraphics(f, fullfile(figDir, ['Raw_vs_' char(primaryMode) '_PeriodHistogram' ext]), ...
            'Resolution', dpi);
        close(f);
    end

    % Figure 3: retention by band.
    RetBand = extended_period_gate_retention(M, "BandName");
    if ~isempty(RetBand)
        f = figure('Color', 'w', 'Position', [100 100 900 650]);
        bar(categorical(RetBand.BandName), RetBand.RetentionPct_EligibleRaw);
        ylabel('Validated retention (% of eligible Raw)', 'FontWeight', 'bold');
        xlabel('Band', 'FontWeight', 'bold');
        title('HSub-supported Raw ultradian candidate retention by band', 'FontWeight', 'bold');
        set(gca, 'FontName', 'Times New Roman', 'Box', 'off', 'TickDir', 'out');
        ylim([0 100]);
        exportgraphics(f, fullfile(figDir, ['Retention_ByBand' ext]), 'Resolution', dpi);
        close(f);
    end

    % Figure 4: retention by photoperiod.
    RetPhoto = extended_period_gate_retention(M, "Photoperiod_h");
    if ~isempty(RetPhoto)
        f = figure('Color', 'w', 'Position', [100 100 900 650]);
        bar(categorical(string(RetPhoto.Photoperiod_h)), RetPhoto.RetentionPct_EligibleRaw);
        ylabel('Validated retention (% of eligible Raw)', 'FontWeight', 'bold');
        xlabel('Photoperiod / light duration (h)', 'FontWeight', 'bold');
        title('HSub-supported Raw ultradian candidate retention by photoperiod', 'FontWeight', 'bold');
        set(gca, 'FontName', 'Times New Roman', 'Box', 'off', 'TickDir', 'out');
        ylim([0 100]);
        exportgraphics(f, fullfile(figDir, ['Retention_ByPhotoperiod' ext]), 'Resolution', dpi);
        close(f);
    end

    % Figure 5: period difference for carried-forward candidates.
    C = M(M.CarryForward & isfinite(M.PeriodDiffPercent), :);
    if ~isempty(C)
        f = figure('Color', 'w', 'Position', [100 100 900 650]);
        scatter(categorical(C.BandName), C.PeriodDiffPercent, 36, 'filled');
        ylabel('Raw vs Selective-HSub period difference (%)', 'FontWeight', 'bold');
        xlabel('Band', 'FontWeight', 'bold');
        title('Period agreement among carried-forward candidates', 'FontWeight', 'bold');
        yline(15, 'k--', '15% tolerance');
        set(gca, 'FontName', 'Times New Roman', 'Box', 'off', 'TickDir', 'out');
        exportgraphics(f, fullfile(figDir, ['CarryForward_PeriodDifference' ext]), 'Resolution', dpi);
        close(f);
    end
end
