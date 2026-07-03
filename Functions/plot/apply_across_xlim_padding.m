function apply_across_xlim_padding(ax, lightHours, varargin)
%APPLY_ACROSS_XLIM_PADDING Offset first/last photoperiod from y-axis spine.
%
%   apply_across_xlim_padding(ax, periodTable.LightDuration_h)
%   apply_across_xlim_padding(ax, lightHours, 'padding', 1)

    if nargin < 2 || isempty(lightHours)
        return;
    end
    if nargin < 1 || isempty(ax)
        ax = gca;
    end

    p = inputParser;
    addParameter(p, 'padding', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});

    hrs = lightHours(:);
    hrs = hrs(isfinite(hrs));
    if isempty(hrs)
        return;
    end

    xlim(ax, [min(hrs) - p.Results.padding, max(hrs) + p.Results.padding]);
end
