function save_user_paths(up)
%SAVE_USER_PATHS Persist user_paths struct to gitignored user_paths.m.
    projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    userFile = fullfile(projectRoot, 'user_paths.m');

    lines = {
        'function up = user_paths()'
        '%USER_PATHS Machine-local output roots (auto-saved).'
        ''
        sprintf('    up.outputRoot = ''%s'';', strrep(up.outputRoot, '''', ''''''))
        sprintf('    up.lastCohort = ''%s'';', strrep(up.lastCohort, '''', ''''''))
        sprintf('    up.lastHSubInput = ''%s'';', strrep(up.lastHSubInput, '''', ''''''))
        sprintf('    up.lastWaveInput = ''%s'';', strrep(up.lastWaveInput, '''', ''''''))
        sprintf('    up.lastAcrossInput = ''%s'';', strrep(up.lastAcrossInput, '''', ''''''))
        'end'
        ''};

    fid = fopen(userFile, 'w');
    if fid < 0
        warning('save_user_paths:WriteFailed', 'Could not write %s', userFile);
        return;
    end
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end
