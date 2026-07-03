function names = prompt_mouse_columns_shared(varNames, defaultIdx)
%PROMPT_MOUSE_COLUMNS_SHARED One-shot mouse column pick for batch runs.

    if ~usejava('desktop')
        names = varNames(defaultIdx);
        return;
    end

    [sel, ok] = listdlg( ...
        'PromptString', 'Select mouse columns (applied to every file in this batch):', ...
        'SelectionMode', 'multiple', ...
        'ListString', varNames, ...
        'InitialValue', defaultIdx, ...
        'ListSize', [320 420]);
    if ~ok || isempty(sel)
        names = {};
    else
        names = varNames(sel);
    end
end
