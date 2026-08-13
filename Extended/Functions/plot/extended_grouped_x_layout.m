function layout = extended_grouped_x_layout(nSubgroups, varargin)
%EXTENDED_GROUPED_X_LAYOUT Centered x-offsets for side-by-side subgroups at each tick.
%
%   layout = extended_grouped_x_layout(nSubgroups)
%   layout = extended_grouped_x_layout(nSubgroups, 'HalfWidth', 0.18, 'Gap', 0.08)
%
%   Returns non-overlapping offsets for nSubgroups (e.g. Female | Male) placed
%   around a categorical x tick. Adjacent box/violin edges are separated by Gap.
%
%   Fields:
%     offsets   — 1×nSubgroups, add to integer group index (1, 2, …)
%     halfWidth — half-width of each subgroup marker/box
%     gap       — horizontal gap between adjacent subgroup edges
%     span      — total width occupied by the group

    p = inputParser;
    addRequired(p, 'nSubgroups', @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'HalfWidth', 0.18, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Gap', 0.08, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, nSubgroups, varargin{:});

    halfW = double(p.Results.HalfWidth);
    gap = double(p.Results.Gap);
    n = max(1, round(double(nSubgroups)));
    step = 2 * halfW + gap;

    offsets = zeros(1, n);
    startC = -((n - 1) * step) / 2;
    for i = 1:n
        offsets(i) = startC + (i - 1) * step;
    end

    layout = struct();
    layout.offsets = offsets;
    layout.halfWidth = halfW;
    layout.gap = gap;
    layout.span = n * 2 * halfW + max(0, n - 1) * gap;
end
