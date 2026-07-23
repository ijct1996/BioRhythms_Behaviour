function paths = setup_paths()
%SETUP_PATHS Add Core v1 + Extended UR project folders and Shared library.
%
%   paths = setup_paths()
%
%   Activates the full BioRhythms_Behaviour tree:
%     Core:     Analysis/, Config/, Functions/, Legacy/
%     Extended: Extended/Analysis/, Config/, Functions/, Legacy/Kent_AG/
%
%   Project paths are prepended so Core/Extended override Shared (e.g. plot_config).
%   setup_extended_paths() remains a thin alias for the same full activation.

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

    %% Core
    addpath(projectRoot, '-begin');
    addpath(paths.analysis, '-begin');
    addpath(paths.config, '-begin');
    if isfolder(paths.legacy)
        addpath(paths.legacy, '-begin');
    end
    if isfolder(paths.functions)
        addpath(genpath(paths.functions), '-begin');
    end

    %% Extended UR (Scripts 4–8 + Kent A–G legacy)
    extendedRoot = fullfile(projectRoot, 'Extended');
    paths.extendedRoot = extendedRoot;
    paths.extendedAnalysis = fullfile(extendedRoot, 'Analysis');
    paths.extendedConfig = fullfile(extendedRoot, 'Config');
    paths.extendedFunctions = fullfile(extendedRoot, 'Functions');
    paths.kentAG = fullfile(extendedRoot, 'Legacy', 'Kent_AG');
    if isfolder(extendedRoot)
        addpath(extendedRoot, '-begin');
        if isfolder(paths.extendedAnalysis)
            addpath(paths.extendedAnalysis, '-begin');
        end
        if isfolder(paths.extendedConfig)
            addpath(paths.extendedConfig, '-begin');
        end
        if isfolder(paths.extendedFunctions)
            addpath(genpath(paths.extendedFunctions), '-begin');
        end
        if isfolder(paths.kentAG)
            addpath(paths.kentAG, '-begin');
        end
        fprintf('BioRhythms Core v1 + Extended UR paths configured.\n');
        fprintf('  Extended root: %s\n', extendedRoot);
    else
        warning('setup_paths:NoExtended', ...
            'Extended folder not found at %s.', extendedRoot);
        fprintf('BioRhythms Core v1 paths configured (Extended missing).\n');
    end
end
