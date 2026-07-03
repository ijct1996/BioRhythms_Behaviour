function [tbl, meta] = read_behav_excel(xlsxPath)
%READ_BEHAV_EXCEL Read behavioural Excel with preserved headers.
    cfg = core_defaults();
    opts = detectImportOptions(xlsxPath, 'FileType', 'spreadsheet');
    opts.VariableNamingRule = 'preserve';
    opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA', '', 'N/A'});
    opts.ImportErrorRule = 'omitrow';
    opts.MissingRule = 'omitrow';
    tbl = readtable(xlsxPath, opts);

    if isempty(tbl) || width(tbl) < 3
        error('read_behav_excel:InsufficientColumns', ...
            'Need Time, >=1 mouse column, and Light duration column.');
    end

    meta = parse_file_metadata(xlsxPath);
    meta.nRows = height(tbl);
    meta.nCols = width(tbl);
    meta.varNames = tbl.Properties.VariableNames;

    [timeIdx, lightIdx, mouseIdx] = preflight_column_mapping(tbl, cfg);
    meta.timeIdx = timeIdx;
    meta.lightIdx = lightIdx;
    meta.mouseIdx = mouseIdx;
    meta.timeName = meta.varNames{timeIdx};
    meta.lightName = meta.varNames{lightIdx};
    meta.mouseNames = meta.varNames(mouseIdx);
    meta.lightDurationHours = light_value_repr(tbl{:, lightIdx});
end
