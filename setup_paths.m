function paths = setup_paths()
%SETUP_PATHS Add project and shared library folders to MATLAB path.
%
%   paths = setup_paths()
%
%   Project paths are prepended so Core v1 overrides Shared (e.g. plot_config).

    thisFile = mfilename('fullpath');
    projectRoot = fileparts(thisFile);
    paths.projectRoot = projectRoot;
    paths.analysis    = fullfile(projectRoot, 'Analysis');
    paths.functions   = fullfile(projectRoot, 'Functions');
    paths.config      = fullfile(projectRoot, 'Config');
    paths.legacy      = fullfile(projectRoot, 'Legacy');
    paths.data        = fullfile(projectRoot, 'Data');
    paths.figures     = fullfile(projectRoot, 'Figures');

    sharedRoot = fullfile(projectRoot, '..', '..', 'Shared');
    if isfolder(sharedRoot)
        paths.sharedRoot = sharedRoot;
        sharedMatlab = fullfile(sharedRoot, 'MATLAB');
        if isfolder(sharedMatlab)
            addpath(genpath(sharedMatlab));
        end
    else
        warning('setup_paths:NoShared', ...
            'Shared library not found at %s.', sharedRoot);
        paths.sharedRoot = '';
    end

    addpath(projectRoot, '-begin');
    addpath(paths.analysis, '-begin');
    addpath(paths.config, '-begin');
    if isfolder(paths.functions)
        addpath(genpath(paths.functions), '-begin');
    end

    fprintf('BioRhythms Core v1 paths configured.\n');
    fprintf('  (Extended UR: run setup_extended_paths from Extended/ — Core paths unchanged.)\n');
end
