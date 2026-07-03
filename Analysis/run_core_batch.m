%% RUN_CORE_BATCH — Script 1 + Script 2 in one unsupervised pass
%
%   All dialogs at start (files, folders, mice, groups).
%   Then harmonic subtraction → wavelet for every selected file.
%
%   setup_paths();
%   run_core_batch

function run_core_batch()
    setup_paths();
    fprintf('BioRhythms Core v1 — Combined batch (HSub + Wavelet)\n');

    cfg = prompt_core_batch_config();
    if cfg.cancelled
        fprintf('Cancelled.\n');
        return;
    end

    nFiles = numel(cfg.files);

    fprintf('=== Phase A: Harmonic subtraction (%d files) ===\n', nFiles);
    for i = 1:nFiles
        inputFile = fullfile(cfg.inPath, cfg.files{i});
        fprintf('\n--- HSub %d/%d: %s ---\n', i, nFiles, cfg.files{i});
        hsub_run_file(inputFile, cfg.hsubRoot, ...
            'Interactive', false, 'MouseNames', cfg.mouseNames, ...
            'PlotMode', cfg.plotMode, 'Groups', cfg.groups, ...
            'SaveIndividualScalograms', cfg.saveIndividualScalograms);
    end

    fprintf('\n=== Phase B: Wavelet (%d files) ===\n', nFiles);
    for i = 1:nFiles
        rawFull = fullfile(cfg.inPath, cfg.files{i});
        fprintf('\n--- Wavelet %d/%d: %s ---\n', i, nFiles, cfg.files{i});
        wavelet_run_file(rawFull, cfg.hsubRoot, cfg.waveRoot, ...
            'PlotMode', cfg.plotMode, 'Interactive', false, ...
            'Groups', cfg.groups);
    end

    fprintf('\nCore batch complete (%d files). Handoff ready for Script 3.\n', nFiles);
end
