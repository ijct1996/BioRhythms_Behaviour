function plotMode = prompt_plot_mode()
%PROMPT_PLOT_MODE Development vs publication figure export (batch setup).
%
%   plotMode = prompt_plot_mode()  % 'development' | 'publication' | '' if cancelled

    up = load_user_paths();
    modes = {'Development (96 DPI PNG)', 'Publication (600 DPI JPEG)'};
    defaultIdx = 1;
    if isfield(up, 'lastPlotMode') && strcmpi(up.lastPlotMode, 'publication')
        defaultIdx = 2;
    end

    [sel, ok] = listdlg('PromptString', 'Plot export mode for this run:', ...
        'SelectionMode', 'single', 'ListString', modes, ...
        'ListSize', [360 90], 'InitialValue', defaultIdx);
    if ~ok
        plotMode = '';
        return;
    end

    if sel == 2
        plotMode = 'publication';
    else
        plotMode = 'development';
    end

    up.lastPlotMode = plotMode;
    save_user_paths(up);
    fprintf('Plot mode: %s\n', plotMode);
end
