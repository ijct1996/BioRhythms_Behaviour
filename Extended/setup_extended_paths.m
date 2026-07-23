function paths = setup_extended_paths()
%SETUP_EXTENDED_PATHS Alias for full Core + Extended path activation.
%
%   paths = setup_extended_paths()
%
%   Prefer setup_paths() from the repo root — it now activates Extended as well.
%   This wrapper remains for older scripts / README snippets that call it.

    thisFile = mfilename('fullpath');
    extendedRoot = fileparts(thisFile);
    projectRoot = fileparts(extendedRoot);

    coreSetup = fullfile(projectRoot, 'setup_paths.m');
    if ~isfile(coreSetup)
        error('setup_extended_paths:NoCoreSetup', ...
            'Core setup_paths.m not found at %s.', coreSetup);
    end

    oldDir = pwd;
    cd(projectRoot);
    cleanupPwd = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    paths = setup_paths();
end
