function groups = group_assignment_dialog(varNames, mouseIdx, meta, varargin)
%GROUP_ASSIGNMENT_DIALOG Legacy-style interactive condition groups.
%
%   Prompts for number of groups, name per group, and mouse columns per group.
%   Matches behav_wavelet_scalogram_v5 / behav_wavelet_powerspectrum_v4 UX.
%
%   groups = group_assignment_dialog(varNames, mouseIdx, meta)
%   groups = group_assignment_dialog(..., 'Interactive', false)  % batch: all mice → one group

    p = inputParser;
    addParameter(p, 'Interactive', true, @islogical);
    parse(p, varargin{:});

    mouseNames = varNames(mouseIdx);

    if ~p.Results.Interactive
        groups = default_single_group(mouseIdx, meta);
        return;
    end

    if ~usejava('desktop')
        warning('group_assignment_dialog:NoDesktop', ...
            'No Java desktop — using single pool for all mouse columns.');
        groups = default_single_group(mouseIdx, meta);
        return;
    end

    prompt = {'Enter number of condition groups:'};
    ansCell = inputdlg(prompt, 'Condition Groups', [1 50]);
    if isempty(ansCell)
        error('group_assignment_dialog:Cancelled', 'Group assignment cancelled.');
    end
    nG = str2double(strtrim(ansCell{1}));
    if isnan(nG) || nG < 1 || mod(nG, 1) ~= 0
        error('group_assignment_dialog:InvalidN', 'Invalid number of groups.');
    end

    groups = struct('name', {}, 'colIdx', {});
    for g = 1:nG
        promptMsg = {sprintf('Enter a name for condition group %d:', g)};
        grpNameCell = inputdlg(promptMsg, 'Condition Group Name', [1 50]);
        if isempty(grpNameCell)
            error('group_assignment_dialog:Cancelled', 'No name for group %d.', g);
        end
        grpName = strtrim(grpNameCell{1});
        if isempty(grpName)
            error('group_assignment_dialog:EmptyName', 'Empty name for group %d.', g);
        end

        [sel, ok] = listdlg( ...
            'PromptString', sprintf('Select activity columns for group "%s":', grpName), ...
            'SelectionMode', 'multiple', ...
            'ListString', mouseNames, ...
            'ListSize', [300 400]);
        if ~ok || isempty(sel)
            error('group_assignment_dialog:NoCols', ...
                'No columns selected for group "%s".', grpName);
        end

        groups(g).name = grpName;
        groups(g).colIdx = mouseIdx(sel);
        groups(g).colNames = mouseNames(sel);
        fprintf('  Group "%s": %d mice\n', grpName, numel(groups(g).colIdx));
    end

    fprintf('Assigned %d condition group(s) for %s\n', nG, meta.fileStem);
end

function groups = default_single_group(mouseIdx, meta)
    if isfield(meta, 'cohort') && strcmp(meta.cohort, 'C57_LP')
        gname = 'C57';
    else
        gname = 'All';
    end
    groups(1).name = gname;
    groups(1).colIdx = mouseIdx(:)';
    groups(1).colNames = meta.varNames(mouseIdx);
    fprintf('Non-interactive mode: single group "%s" (%d mice).\n', gname, numel(mouseIdx));
end
