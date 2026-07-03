function info = parse_mouse_metadata(mouseNames)
%PARSE_MOUSE_METADATA Extract sex/genotype hints from column headers.
%
%   Recognises suffixes: -M, -F, _M, _F, Male, Female (case-insensitive).

    n = numel(mouseNames);
    info = struct('name', mouseNames(:), 'sex', repmat({''}, n, 1), ...
        'genotypeHint', repmat({''}, n, 1));

    for i = 1:n
        nm = char(mouseNames{i});
        info(i).sex = infer_sex(nm);
        if contains(upper(nm), 'NR2B')
            info(i).genotypeHint = 'NR2B';
        elseif contains(upper(nm), 'C57') || contains(upper(nm), 'WT')
            info(i).genotypeHint = 'C57';
        end
    end
end

function sex = infer_sex(nm)
    sex = '';
    nmLower = lower(nm);
    isMale = ~isempty(regexp(nm, '(?:^|[-_])M(?:$|[-_])', 'once', 'ignorecase')) || ...
        (contains(nmLower, 'male') && ~contains(nmLower, 'female'));
    isFemale = ~isempty(regexp(nm, '(?:^|[-_])F(?:$|[-_])', 'once', 'ignorecase')) || ...
        contains(nmLower, 'female');
    if isMale
        sex = 'M';
    elseif isFemale
        sex = 'F';
    end
end
