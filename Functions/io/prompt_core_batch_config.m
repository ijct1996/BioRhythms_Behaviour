function cfg = prompt_core_batch_config()
%PROMPT_CORE_BATCH_CONFIG One-shot setup for Script 1 + Script 2 batch run.

    up = load_user_paths();
    startDir = up.lastHSubInput;
    if isempty(startDir) || ~isfolder(startDir)
        startDir = up.lastWaveInput;
    end
    if isempty(startDir) || ~isfolder(startDir)
        startDir = pwd;
    end

    cfg = struct();
    cfg.cancelled = false;

    fprintf('=== Core v1 batch setup (all prompts now; then unsupervised) ===\n');

    [cfg.files, cfg.inPath, cfg.inputMode] = prompt_input_files(startDir, ...
        'Select Pipeline_Input xlsx file(s)');
    if isempty(cfg.files)
        cfg.cancelled = true;
        return;
    end
    fprintf('Selected %d file(s) [%s mode]\n', numel(cfg.files), cfg.inputMode);

    cfg.hsubRoot = uigetdir(cfg.inPath, ...
        'Select/create 01_HarmonicSubtraction folder');
    if isequal(cfg.hsubRoot, 0)
        cfg.cancelled = true;
        return;
    end
    ensure_dir(cfg.hsubRoot);

    firstFile = fullfile(cfg.inPath, cfg.files{1});
    [~, meta] = read_behav_excel(firstFile);

    cfg.mouseNames = prompt_mouse_columns_shared(meta.varNames, meta.mouseIdx);
    if isempty(cfg.mouseNames)
        cfg.cancelled = true;
        return;
    end

    defaultWave = fullfile(fileparts(cfg.hsubRoot), '02_Wavelet');
    if isfolder(defaultWave)
        waveStart = defaultWave;
    else
        waveStart = fileparts(cfg.hsubRoot);
    end
    cfg.waveRoot = uigetdir(waveStart, ...
        'Select/create 02_Wavelet folder');
    if isequal(cfg.waveRoot, 0)
        cfg.cancelled = true;
        return;
    end
    ensure_dir(cfg.waveRoot);

    groups = group_assignment_dialog(meta.varNames, meta.mouseIdx, meta, ...
        'Interactive', true);
    cfg.groups = groups_to_template(groups, meta.varNames);

    cfg.plotMode = prompt_plot_mode();
    if isempty(cfg.plotMode)
        cfg.cancelled = true;
        return;
    end

    cfg.saveIndividualScalograms = prompt_hsub_individual_scalograms();

    fprintf('\n--- Batch config summary ---\n');
    fprintf('  Files     : %d\n', numel(cfg.files));
    fprintf('  Mice      : %d columns\n', numel(cfg.mouseNames));
    fprintf('  Groups    : %d\n', numel(cfg.groups));
    fprintf('  Plot mode : %s\n', cfg.plotMode);
    fprintf('  HSub indiv: %s\n', onoff(cfg.saveIndividualScalograms));
    fprintf('  HSub out  : %s\n', cfg.hsubRoot);
    fprintf('  Wave out  : %s\n', cfg.waveRoot);
    fprintf('  Handoff   : %s\n', fullfile(fileparts(cfg.waveRoot), 'Handoff'));
    fprintf('--- Starting unsupervised run ---\n\n');

    up.lastHSubInput = cfg.inPath;
    up.lastWaveInput = cfg.inPath;
    save_user_paths(up);
end

function s = onoff(tf)
    if tf, s = 'yes (150 DPI PNG)'; else, s = 'no'; end
end
