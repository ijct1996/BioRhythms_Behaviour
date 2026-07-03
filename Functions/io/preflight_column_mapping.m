function [timeIdx, lightIdx, mouseIdx] = preflight_column_mapping(tbl, cfg)
%PREFLIGHT_COLUMN_MAPPING Infer Time, Light, and mouse column indices.

    if nargin < 2
        cfg = core_defaults();
    end

    names = tbl.Properties.VariableNames;
    timeIdx = find_first_name(names, cfg.io.timeColumnCandidates, 1);
    lightIdx = find_first_name(names, cfg.io.lightColumnCandidates, width(tbl));
    mouseIdx = setdiff(1:width(tbl), [timeIdx, lightIdx], 'stable');

    if isempty(mouseIdx)
        error('preflight_column_mapping:NoMice', 'No mouse columns found.');
    end
end

function idx = find_first_name(names, candidates, fallback)
    idx = [];
    for k = 1:numel(candidates)
        hit = find(strcmpi(strtrim(names), strtrim(candidates{k})), 1);
        if ~isempty(hit)
            idx = hit;
            return;
        end
    end
    idx = fallback;
end
