function update_handoff_index(handoffDir, row)
%UPDATE_HANDOFF_INDEX Append row to CoreHandoff_Index.xlsx.

    cfg = core_defaults();
    indexPath = fullfile(handoffDir, cfg.handoff.indexName);
    ensure_dir(handoffDir);

    if isfile(indexPath)
        idxTbl = readtable(indexPath);
        idxTbl = [idxTbl; row]; %#ok<AGROW>
    else
        idxTbl = row;
    end

    writetable(idxTbl, indexPath);
end
