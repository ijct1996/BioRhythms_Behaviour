function saveIndividual = prompt_hsub_individual_scalograms()
%PROMPT_HSUB_INDIVIDUAL_SCALOGRAMS Optional per-mouse HSub scalograms (150 DPI PNG).

    choices = { ...
        'No — group-average HSub scalograms only (recommended)', ...
        'Yes — also save individual HSub scalograms (150 DPI PNG)'};
    [sel, ok] = listdlg('PromptString', ...
        'Individual HSub scalogram export (QC / supplementary):', ...
        'SelectionMode', 'single', 'ListString', choices, 'ListSize', [420 90]);
    saveIndividual = ok && sel == 2;
    if saveIndividual
        fprintf('Individual HSub scalograms: ON (150 DPI PNG)\n');
    else
        fprintf('Individual HSub scalograms: OFF\n');
    end
end
