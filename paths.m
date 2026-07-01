function paths = setup_paths()
%SETUP_PATHS Add project and shared library folders to MATLAB path.
%
%   paths = setup_paths()
%
%   Resolves paths relative to this file — no hard-coded absolute paths.
%   Call at the start of each session or from startup.m.
%
%   Outputs:
%       paths — struct with resolved directory paths
%
%   See also: Shared/Research_Knowledge.md

    %% Project root
    thisFile = mfilename('fullpath');
    projectRoot = fileparts(thisFile);
    paths.projectRoot = projectRoot;

    %% Project subfolders
    paths.analysis   = fullfile(projectRoot, 'Analysis');
    paths.functions  = fullfile(projectRoot, 'Functions');
    paths.figures    = fullfile(projectRoot, 'Figures');
    paths.data       = fullfile(projectRoot, 'Data');

    addpath(paths.analysis);
    addpath(paths.functions);

    %% Shared library (Option B: sibling clone at Research/Shared/)
  sharedRoot = fullfile(projectRoot, '..', '..', 'Shared');
    if isfolder(sharedRoot)
        paths.sharedRoot = sharedRoot;
        addpath(genpath(fullfile(sharedRoot, 'MATLAB')));
    else
        warning('setup_paths:NoShared', ...
            'Shared library not found at %s. Clone research-shared or add submodule.', sharedRoot);
        paths.sharedRoot = '';
    end

    fprintf('Paths configured. Project: %s\n', projectRoot);
end
