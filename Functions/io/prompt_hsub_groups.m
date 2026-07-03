function groups = prompt_hsub_groups(varNames, mouseIdx, meta)
%PROMPT_HSUB_GROUPS Condition groups or cohort-wide average for HSub scalograms.

    choices = { ...
        'Assign condition groups (Male / Female / …)', ...
        'Single cohort-average group (All)'};
    [sel, ok] = listdlg('PromptString', ...
        'Group assignment for HSub scalogram averages:', ...
        'SelectionMode', 'single', 'ListString', choices, 'ListSize', [380 90]);
    if ~ok
        groups = [];
        return;
    end

    if sel == 1
        groups = group_assignment_dialog(varNames, mouseIdx, meta, 'Interactive', true);
    else
        groups = group_assignment_dialog(varNames, mouseIdx, meta, 'Interactive', false);
    end
    groups = groups_to_template(groups, varNames);
end
