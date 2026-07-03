function outPaths = across_plot_period_distributions(periodTable, figFolder, baseName, theme, varargin)
%ACROSS_PLOT_PERIOD_DISTRIBUTIONS Per-sex box + jittered points at each photoperiod.
%
%   One figure per condition group (e.g. F, M). Box = median/IQR; points = mice.
%   Recommended over combined scatter when rank-1 UR periods overlap vertically.

    if nargin < 4 || isempty(theme)
        theme = plot_config('development');
    end

    p = inputParser;
    addParameter(p, 'titlePrefix', 'UR period distribution', @ischar);
    addParameter(p, 'ylabelStr', 'Validated UR period (h) — rank 1', @ischar);
    addParameter(p, 'ylimRange', [0 18], @isnumeric);
    addParameter(p, 'jitterWidth', 0.35, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});

    outPaths = {};
    if isempty(periodTable) || height(periodTable) == 0
        warning('across_plot_period_distributions:Empty', 'No period data.');
        return;
    end

    top1 = periodTable(periodTable.PeakRank == 1, :);
    if isempty(top1), return; end

    distFolder = fullfile(figFolder, 'Period_Distributions');
    ensure_dir(distFolder);
    ext = theme.exportFormat;
    groups = unique(top1.Group);

    for g = 1:numel(groups)
        grpTag = groups(g);
        sub = top1(top1.Group == grpTag, :);
        if isempty(sub), continue; end

        safeTag = sanitise_filename(char(grpTag));
        outFile = fullfile(distFolder, sprintf('%s_%s.%s', baseName, safeTag, ext));
        across_plot_period_distribution_one(sub, outFile, grpTag, theme, p.Results);
        outPaths{end + 1} = outFile; %#ok<AGROW>
        fprintf('  Period distribution: %s → %s\n', collaborator_group_label(grpTag), outFile);
    end
end

function across_plot_period_distribution_one(sub, outPath, grpTag, theme, opts)

    photos = sort(unique(sub.LightDuration_h));
    photos = photos(isfinite(photos));
    if isempty(photos), return; end

    st = pick_collaborator_group_style(grpTag, theme.palette);
    sexLabel = collaborator_group_label(grpTag);

    xAll = [];
    yAll = [];
    for k = 1:numel(photos)
        yk = sub.PeakPeriod_hr(sub.LightDuration_h == photos(k));
        yk = yk(isfinite(yk));
        xAll = [xAll; repmat(photos(k), numel(yk), 1)]; %#ok<AGROW>
        yAll = [yAll; yk(:)]; %#ok<AGROW>
    end
    if isempty(yAll), return; end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 960 520]);
    ax = axes(fig);
    hold(ax, 'on');

    boxchart(ax, xAll, yAll, ...
        'BoxFaceColor', st.lineColor, ...
        'BoxFaceAlpha', 1, ...
        'BoxEdgeColor', [0 0 0], ...
        'LineWidth', 1, ...
        'MarkerStyle', 'none', ...
        'BoxWidth', 0.55);

    rng(42);
    ptSize = 36;
    for k = 1:numel(photos)
        L = photos(k);
        yk = sub.PeakPeriod_hr(sub.LightDuration_h == L);
        yk = yk(isfinite(yk));
        n = numel(yk);
        if n == 0, continue; end
        jitter = (rand(n, 1) - 0.5) * 2 * opts.jitterWidth;
        scatter(ax, L + jitter, yk, ptSize, ...
            'Marker', st.marker, ...
            'MarkerFaceColor', st.markerFaceColor, ...
            'MarkerEdgeColor', st.markerEdgeColor, ...
            'LineWidth', 1);
    end
    hold(ax, 'off');

    xticks(ax, photos);
    xlabel(ax, 'Light duration (h)');
    ylabel(ax, opts.ylabelStr);
    title(ax, sprintf('%s — %s', opts.titlePrefix, sexLabel), 'Interpreter', 'none');
    ylim(ax, opts.ylimRange);
    apply_across_xlim_padding(ax, photos);
    apply_across_figure_theme(fig, theme);
    export_figure(fig, outPath, theme);
    close(fig);
end

function s = sanitise_filename(strIn)
    s = char(regexprep(string(strIn), '[^\w\-]', '_'));
end
