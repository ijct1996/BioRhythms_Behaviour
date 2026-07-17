function outDir = extended_plot_transition_gradient(resync, cfg)
%EXTENDED_PLOT_TRANSITION_GRADIENT Photoperiod × transition effect-size summary.
%
%   Co-primary endpoints: DeltaR (phase coherence) and ridge-power change.
%   Primary bands: UR_1_3 and UR_3_6 (cfg.bands.primaryUR).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    outDir = '';
    if isempty(resync) || ~isstruct(resync)
        return;
    end

    if isfield(resync, 'outRoot')
        outRoot = resync.outRoot;
    elseif isfield(resync, 'outXLSX')
        outRoot = fileparts(resync.outXLSX);
    else
        return;
    end

    outDir = fullfile(outRoot, 'TransitionEffect_vs_Photoperiod');
    extended_period_gate_ensure_dir(outDir);

    primaryBands = string(cfg.bands.primaryUR);
    transOrder = ["DL", "LD"];
    dpi = cfg.plot.saveDpi;
    ext = cfg.plot.figExt;

    % --- Build summary table from CandidateDeltaR and ridge power stats ---
    gradRows = {};
    hdr = {'Photoperiod_h', 'BandName', 'TransitionType', 'Metric', ...
        'MeanEffect', 'MedianEffect', 'SD_Effect', 'N_Candidates', 'N_Significant_BH'};

    if isfield(resync, 'CandidateDeltaR') && ~isempty(resync.CandidateDeltaR)
        CD = resync.CandidateDeltaR;
        CD = CD(ismember(string(CD.TransitionType), transOrder) & ...
            ismember(string(CD.BandName), primaryBands), :);
        groups = unique(CD(:, {'Photoperiod_h', 'BandName', 'TransitionType'}), 'rows', 'stable');
        for g = 1:height(groups)
            idx = CD.Photoperiod_h == groups.Photoperiod_h(g) & ...
                string(CD.BandName) == string(groups.BandName(g)) & ...
                string(CD.TransitionType) == string(groups.TransitionType(g));
            x = CD.DeltaR_PostMinusPre(idx);
            x = x(isfinite(x));
            gradRows(end+1, :) = {groups.Photoperiod_h(g), string(groups.BandName(g)), ...
                string(groups.TransitionType(g)), "DeltaR", ...
                mean(x, 'omitnan'), median(x, 'omitnan'), std(x, 'omitnan'), numel(x), NaN}; %#ok<AGROW>
        end
    end

    if isfield(resync, 'Resync_RidgePowerStats_BH_FDR') && ~isempty(resync.Resync_RidgePowerStats_BH_FDR)
        RP = resync.Resync_RidgePowerStats_BH_FDR;
        RP = RP(ismember(string(RP.TransitionType), transOrder) & ...
            ismember(string(RP.BandName), primaryBands), :);
        if ismember('MeanDifference_PostMinusPre', RP.Properties.VariableNames)
            metricCol = 'MeanDifference_PostMinusPre';
        else
            metricCol = '';
        end
        if ~isempty(metricCol)
            groups = unique(RP(:, {'Photoperiod_h', 'BandName', 'TransitionType'}), 'rows', 'stable');
            for g = 1:height(groups)
                idx = RP.Photoperiod_h == groups.Photoperiod_h(g) & ...
                    string(RP.BandName) == string(groups.BandName(g)) & ...
                    string(RP.TransitionType) == string(groups.TransitionType(g));
                x = RP.(metricCol)(idx);
                x = x(isfinite(x));
                nSig = 0;
                if ismember('Significant_BH', RP.Properties.VariableNames)
                    nSig = sum(logical(RP.Significant_BH(idx)));
                end
                gradRows(end+1, :) = {groups.Photoperiod_h(g), string(groups.BandName(g)), ...
                    string(groups.TransitionType(g)), "RidgePower_PostMinusPre", ...
                    mean(x, 'omitnan'), median(x, 'omitnan'), std(x, 'omitnan'), numel(x), nSig}; %#ok<AGROW>
            end
        end
    end

    if isempty(gradRows)
        writecell({'No gradient rows.'}, fullfile(outDir, 'TransitionEffect_vs_Photoperiod.xlsx'), ...
            'Sheet', 'Summary');
        return;
    end

    Grad = cell2table(gradRows, 'VariableNames', hdr);
    outXlsx = fullfile(outDir, 'TransitionEffect_vs_Photoperiod.xlsx');
    writetable(Grad, outXlsx, 'Sheet', 'Summary');

    if ~cfg.plot.generateFigures
        return;
    end

    photos = sort(unique(Grad.Photoperiod_h));
    photoLabels = strings(size(photos));
    for i = 1:numel(photos)
        if photos(i) >= 24
            photoLabels(i) = "LL";
        else
            photoLabels(i) = "L" + string(photos(i));
        end
    end

    metrics = unique(string(Grad.Metric), 'stable');
    for mi = 1:numel(metrics)
        metric = metrics(mi);
        Gm = Grad(string(Grad.Metric) == metric, :);
        for bi = 1:numel(primaryBands)
            band = primaryBands(bi);
            Gb = Gm(string(Gm.BandName) == band, :);
            if isempty(Gb), continue; end

            valsDL = nan(numel(photos), 1);
            valsLD = nan(numel(photos), 1);
            for p = 1:numel(photos)
                rowDL = Gb(Gb.Photoperiod_h == photos(p) & string(Gb.TransitionType) == "DL", :);
                rowLD = Gb(Gb.Photoperiod_h == photos(p) & string(Gb.TransitionType) == "LD", :);
                if ~isempty(rowDL), valsDL(p) = rowDL.MeanEffect(1); end
                if ~isempty(rowLD), valsLD(p) = rowLD.MeanEffect(1); end
            end

            f = figure('Color', 'w', 'Position', [100 100 1000 520]);
            hold on;
            plot(1:numel(photos), valsDL, '-o', 'LineWidth', 1.5, 'DisplayName', 'DL (lights-on)');
            plot(1:numel(photos), valsLD, '-s', 'LineWidth', 1.5, 'DisplayName', 'LD (lights-off)');
            yline(0, 'k:', 'LineWidth', 1);
            xticks(1:numel(photos));
            xticklabels(cellstr(photoLabels));
            xlabel('Photoperiod');
            ylabel(char(metric));
            title(sprintf('%s | %s across photoperiod', char(metric), char(band)), 'Interpreter', 'none');
            legend('Location', 'best');
            box off;
            set(gca, 'TickDir', 'out');
            safeName = regexprep(char(band), '[^\w]', '_');
            exportgraphics(f, fullfile(outDir, sprintf('Gradient_%s_%s%s', char(metric), safeName, ext)), ...
                'Resolution', dpi);
            close(f);
        end
    end

    fprintf('Transition gradient summary: %s\n', outXlsx);
end
