function groups = resolve_groups_by_names(varNames, groupsTemplate)
%RESOLVE_GROUPS_BY_NAMES Apply saved group template (by column name) to a new file.

    groups = struct('name', {}, 'colIdx', {}, 'colNames', {});
    for g = 1:numel(groupsTemplate)
        if isfield(groupsTemplate(g), 'colNames') && ~isempty(groupsTemplate(g).colNames)
            names = groupsTemplate(g).colNames;
            if isstring(names), names = cellstr(names); end
        else
            names = varNames(groupsTemplate(g).colIdx);
        end
        groups(g).name = groupsTemplate(g).name;
        groups(g).colNames = names;
        groups(g).colIdx = resolve_column_indices(varNames, names);
    end
end
