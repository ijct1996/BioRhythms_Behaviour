function meta = parse_collaborator_group(grp)
%PARSE_COLLABORATOR_GROUP Sex and genotype from condition group label.
%
%   meta.sex      — 'F', 'M', or ''
%   meta.genotype — 'Ctrl', 'KO', or ''
%
%   Examples: "F - Ctrl", "F - KO", "M - Ctrl", "Male", "Female"

    g = lower(strtrim(char(grp)));
    meta = struct('sex', '', 'genotype', '', 'raw', g);

    if contains(g, 'female') || ~isempty(regexp(g, '^f(\s|$|-|–)', 'once'))
        meta.sex = 'F';
    elseif ~isempty(regexp(g, '^m(\s|$|-|–)', 'once')) || ...
            (contains(g, 'male') && ~contains(g, 'female'))
        meta.sex = 'M';
    end

    if contains(g, 'ko') || contains(g, 'knockout') || contains(g, '-/-')
        meta.genotype = 'KO';
    elseif contains(g, 'ctrl') || contains(g, 'control') || ...
            contains(g, 'wt') || contains(g, 'wild') || contains(g, '+/+')
        meta.genotype = 'Ctrl';
    end
end
