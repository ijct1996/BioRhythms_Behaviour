function [files, folder, mode] = prompt_input_files(startDir, dlgTitle)
%PROMPT_INPUT_FILES Legacy-style single vs multiple Excel file picker.
%
%   [files, folder, mode] = prompt_input_files()
%   mode — 'single' | 'multiple'

    if nargin < 1 || isempty(startDir) || ~isfolder(startDir)
        startDir = pwd;
    end
    if nargin < 2 || isempty(dlgTitle)
        dlgTitle = 'Select input Excel file(s)';
    end

    mode = '';
    files = {};
    folder = '';

    if usejava('desktop')
        choice = questdlg('Process one file or multiple files?', ...
            'Input mode', 'Single file', 'Multiple files', 'Multiple files');
        if isempty(choice)
            return;
        end
        multi = strcmp(choice, 'Multiple files');
        mode = lower(strrep(choice, ' file', ''));
    else
        multi = true;
        mode = 'multiple';
    end

    if multi
        multiArg = 'on';
    else
        multiArg = 'off';
    end

    [f, folder] = uigetfile({'*.xlsx', 'Excel (*.xlsx)'}, dlgTitle, startDir, ...
        'MultiSelect', multiArg);
    if isequal(f, 0)
        files = {};
        folder = '';
        return;
    end

    if ischar(f)
        files = {f};
    else
        files = cellstr(f);
    end
    files = sort(files(:))';
end
