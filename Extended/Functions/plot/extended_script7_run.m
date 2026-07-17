function out = extended_script7_run(extHandoffDir, cfg)
%EXTENDED_SCRIPT7_RUN Phase events + publication profiles (Kent F + G merged).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    fprintf('\n=== Extended Script 7: Phase events + profiles (F+G) ===\n');

    if ~cfg.plot.generateLegacyFigures
        fprintf(['Skipped legacy Kent F/G in development mode (generateLegacyFigures=false).\n' ...
            '  Enable via: cfg = extended_defaults(); cfg.plot.generateLegacyFigures=true;\n' ...
            '  Or request publication pass later.\n']);
        out = struct('skipped', true, 'reason', 'generateLegacyFigures=false');
        return;
    end

    fprintf('Warning: Kent F/G still output 600 DPI JPEG until modularised.\n');

    kentDir = extended_kent_ag_dir_();
    scripts = { ...
        'F_Plot_ValidatedUR_PhaseEventHistograms_v1.m', ...
        'G_Plot_SelectedValidatedUR_PublicationProfiles_v1.m'};

    addpath(kentDir, '-begin');
    oldDir = pwd;
    cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(kentDir);

    for i = 1:numel(scripts)
        scriptPath = fullfile(kentDir, scripts{i});
        if ~isfile(scriptPath)
            error('extended_script7_run:MissingScript', 'Not found: %s', scriptPath);
        end
        fprintf('Running %s ...\n', scripts{i});
        run(scripts{i});
    end

    out = struct('scripts', {scripts}, 'kentDir', kentDir);
    fprintf('\nExtended Script 7 complete (legacy Kent F+G).\n');
end

function kentDir = extended_kent_ag_dir_()
    thisFile = mfilename('fullpath');
    plotDir = fileparts(thisFile);
    funcDir = fileparts(plotDir);
    extendedRoot = fileparts(funcDir);
    kentDir = fullfile(extendedRoot, 'Legacy', 'Kent_AG');
end
