function template = groups_to_template(groups, varNames)
%GROUPS_TO_TEMPLATE Store column names for batch reuse across photoperiod files.

    template = groups;
    for g = 1:numel(template)
        if ~isfield(template(g), 'colNames') || isempty(template(g).colNames)
            template(g).colNames = varNames(template(g).colIdx);
        end
    end
end
