function paths = setup_extended_paths()
%SETUP_EXTENDED_PATHS Add Extended UR folders and Core project paths.
%
%   paths = setup_extended_paths()
%
%   Prefer this entry for Extended work. Core setup_paths() behaviour is
%   unchanged aside from a one-line note about Extended/.
%
%   See also: Extended/README.md, setup_paths

    thisFile = mfilename('fullpath');
    extendedRoot = fileparts(thisFile);
    projectRoot = fileparts(extendedRoot);

    %% Core paths first (Shared + Analysis/Functions/Config)
    coreSetup = fullfile(projectRoot, 'setup_paths.m');
    if isfile(coreSetup)
        oldDir = pwd;
        cd(projectRoot);
        cleanupPwd = onCleanup(@() cd(oldDir));
        corePaths = setup_paths();
        clear cleanupPwd;
    else
        warning('setup_extended_paths:NoCoreSetup', ...
            'Core setup_paths.m not found at %s.', coreSetup);
        corePaths = struct('projectRoot', projectRoot);
    end

    paths = corePaths;
    paths.extendedRoot = extendedRoot;
    paths.extendedAnalysis = fullfile(extendedRoot, 'Analysis');
    paths.extendedConfig = fullfile(extendedRoot, 'Config');
    paths.extendedFunctions = fullfile(extendedRoot, 'Functions');
    paths.kentAG = fullfile(extendedRoot, 'Legacy', 'Kent_AG');

    addpath(extendedRoot, '-begin');
    addpath(paths.extendedAnalysis, '-begin');
    addpath(paths.extendedConfig, '-begin');
    if isfolder(paths.extendedFunctions)
        addpath(genpath(paths.extendedFunctions), '-begin');
    end
    if isfolder(paths.kentAG)
        addpath(paths.kentAG, '-begin');
    end

    fprintf('Extended UR paths configured (%s).\n', extendedRoot);
end
