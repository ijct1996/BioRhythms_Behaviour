function out = wavelet_vertcat_tables(cells)
%WAVELET_VERTCAT_TABLES Vertically concat tables with differing columns (skip empty).

    if isempty(cells)
        out = table();
        return;
    end

    keep = false(numel(cells), 1);
    for i = 1:numel(cells)
        keep(i) = istable(cells{i}) && height(cells{i}) > 0;
    end
    cells = cells(keep);
    if isempty(cells)
        out = table();
        return;
    end

    allVars = cells{1}.Properties.VariableNames;
    for i = 2:numel(cells)
        allVars = union(allVars, cells{i}.Properties.VariableNames, 'stable');
    end

    for i = 1:numel(cells)
        T = cells{i};
        missingVars = setdiff(allVars, T.Properties.VariableNames, 'stable');
        for v = 1:numel(missingVars)
            T.(missingVars{v}) = default_missing_column(T, missingVars{v});
        end
        cells{i} = T(:, allVars);
    end

    out = vertcat(cells{:});
end

function col = default_missing_column(T, varName)
    n = height(T);
    if n == 0
        col = [];
        return;
    end
    if endsWith(varName, {'ID', 'Group'}, 'IgnoreCase', true) || contains(varName, 'Signal')
        col = strings(n, 1);
        col(:) = missing;
    else
        col = nan(n, 1);
    end
end
