%% RUN_HARMONIC_SUBTRACTION — Core v1 Script 1
% All dialogs at start; batch runs unsupervised.
%
%   setup_paths();
%   run_harmonic_subtraction

function run_harmonic_subtraction()
    setup_paths();
    fprintf('BioRhythms Core v1 — Harmonic subtraction\n');

    cfg = prompt_hsub_batch_config();
    if cfg.cancelled
        fprintf('Cancelled.\n');
        return;
    end

    nFiles = numel(cfg.files);
    fprintf('\n--- Running HSub on %d file(s) (no further prompts) ---\n', nFiles);

    for i = 1:nFiles
        inputFile = fullfile(cfg.inPath, cfg.files{i});
        fprintf('\n=== HSub %d/%d: %s ===\n', i, nFiles, cfg.files{i});
        hsub_run_file(inputFile, cfg.outRoot, ...
            'Interactive', false, 'MouseNames', cfg.mouseNames, ...
            'PlotMode', cfg.plotMode, 'Groups', cfg.groups, ...
            'SaveIndividualScalograms', cfg.saveIndividualScalograms);
    end

    fprintf('\nAll harmonic subtraction runs complete (%d files).\n', nFiles);
end
