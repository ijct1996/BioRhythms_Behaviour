function label = collaborator_group_label(grp)
%COLLABORATOR_GROUP_LABEL Display name for sex/genotype group tags.

    g = lower(strtrim(char(grp)));
    if is_collaborator_female_group(g)
        label = 'Female';
    elseif is_collaborator_male_group(g)
        label = 'Male';
    else
        label = char(grp);
    end
end

function tf = is_collaborator_female_group(g)
    tf = any(strcmp(g, {'f', 'female', 'fem', 'females'})) || contains(g, 'female');
end

function tf = is_collaborator_male_group(g)
    tf = any(strcmp(g, {'m', 'male', 'males'})) ...
        || (contains(g, 'male') && ~contains(g, 'female'));
end
