function idx = resolve_column_indices(varNames, colNames)
%RESOLVE_COLUMN_INDICES Map column header names to table indices (batch reuse).

    if isstring(colNames)
        colNames = cellstr(colNames);
    end
    idx = zeros(1, numel(colNames));
    for i = 1:numel(colNames)
        hit = find(strcmp(varNames, colNames{i}), 1);
        if isempty(hit)
            error('resolve_column_indices:MissingColumn', ...
                'Column "%s" not found in this file.', colNames{i});
        end
        idx(i) = hit;
    end
end
