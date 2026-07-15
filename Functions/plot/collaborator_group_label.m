function label = collaborator_group_label(grp)
%COLLABORATOR_GROUP_LABEL Human-readable label for legend / filenames.

    meta = parse_collaborator_group(grp);
    if strcmp(meta.sex, 'F') && strcmp(meta.genotype, 'Ctrl')
        label = 'Female Control';
    elseif strcmp(meta.sex, 'F') && strcmp(meta.genotype, 'KO')
        label = 'Female KO';
    elseif strcmp(meta.sex, 'M') && strcmp(meta.genotype, 'Ctrl')
        label = 'Male Control';
    elseif strcmp(meta.sex, 'M') && strcmp(meta.genotype, 'KO')
        label = 'Male KO';
    elseif strcmp(meta.sex, 'F')
        label = 'Female';
    elseif strcmp(meta.sex, 'M')
        label = 'Male';
    else
        label = char(grp);
    end
end
