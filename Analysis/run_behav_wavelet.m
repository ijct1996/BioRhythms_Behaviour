%% RUN_BEHAV_WAVELET — Core v1 Script 2
% All dialogs at start; batch runs unsupervised.
%
%   setup_paths();
%   run_behav_wavelet

function run_behav_wavelet()
    setup_paths();
    fprintf('BioRhythms Core v1 — Wavelet\n');

    cfg = prompt_wavelet_batch_config();
    if cfg.cancelled
        fprintf('Cancelled.\n');
        return;
    end

    nFiles = numel(cfg.files);
    fprintf('\n--- Running wavelet on %d file(s) (no further prompts) ---\n', nFiles);

    for i = 1:nFiles
        rawFull = fullfile(cfg.rawPath, cfg.files{i});
        fprintf('\n=== Wavelet %d/%d: %s ===\n', i, nFiles, cfg.files{i});
        wavelet_run_file(rawFull, cfg.hsubRoot, cfg.outRoot, ...
            'PlotMode', cfg.plotMode, 'Interactive', false, ...
            'Groups', cfg.groups);
    end

    fprintf('\nAll wavelet runs complete (%d files).\n', nFiles);
end
