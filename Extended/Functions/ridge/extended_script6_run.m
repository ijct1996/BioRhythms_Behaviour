function out = extended_script6_run(extHandoffDir, cfg)
%EXTENDED_SCRIPT6_RUN Across-photoperiod LME / FDR on validated raw UR (Kent E).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    fprintf('\n=== Extended Script 6: Across-photoperiod LME / FDR ===\n');
    fprintf('Note: legacy Kent E may still emit 600 DPI JPEG figures until modularised.\n');
    fprintf('      LME/FDR tables are the primary development deliverable.\n');

    kentDir = extended_kent_ag_dir_();
    scriptName = 'E_AcrossPhotoperiod_Analysis_v8_1_ValidatedRawUR_FDR_SexStable_PCAFix.m';
    scriptPath = fullfile(kentDir, scriptName);
    if ~isfile(scriptPath)
        error('extended_script6_run:MissingScript', 'Not found: %s', scriptPath);
    end

    if nargin >= 1 && ~isempty(extHandoffDir)
        fprintf('Extended handoff context: %s\n', extHandoffDir);
    end

    addpath(kentDir, '-begin');
    oldDir = pwd;
    cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(kentDir);
    run(scriptName);

    out = struct('script', scriptName, 'kentDir', kentDir);
    fprintf('\nExtended Script 6 complete (legacy Kent E).\n');
end

function kentDir = extended_kent_ag_dir_()
    thisFile = mfilename('fullpath');
    ridgeDir = fileparts(thisFile);
    funcDir = fileparts(ridgeDir);
    extendedRoot = fileparts(funcDir);
    kentDir = fullfile(extendedRoot, 'Legacy', 'Kent_AG');
end
