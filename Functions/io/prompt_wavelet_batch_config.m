function cfg = prompt_wavelet_batch_config()
%PROMPT_WAVELET_BATCH_CONFIG Collect all Script 2 settings upfront (unsupervised batch).

    up = load_user_paths();
    startDir = up.lastWaveInput;
    if isempty(startDir) || ~isfolder(startDir)
        startDir = pwd;
    end

    cfg = struct();
    cfg.version = core_defaults().version;

    [cfg.files, cfg.rawPath, cfg.inputMode] = prompt_input_files(startDir, ...
        'Select RAW Pipeline_Input xlsx file(s)');
    if isempty(cfg.files)
        cfg.cancelled = true;
        return;
    end
    cfg.cancelled = false;
    fprintf('Wavelet batch: %s mode, %d file(s)\n', cfg.inputMode, numel(cfg.files));

    cfg.hsubRoot = uigetdir(cfg.rawPath, ...
        'Select 01_HarmonicSubtraction folder for this cohort');
    if isequal(cfg.hsubRoot, 0)
        cfg.cancelled = true;
        return;
    end

    cfg.outRoot = uigetdir(cfg.hsubRoot, ...
        'Select/create 02_Wavelet output folder');
    if isequal(cfg.outRoot, 0)
        cfg.cancelled = true;
        return;
    end
    ensure_dir(cfg.outRoot);

    firstFile = fullfile(cfg.rawPath, cfg.files{1});
    [~, meta] = read_behav_excel(firstFile);
    groups = group_assignment_dialog(meta.varNames, meta.mouseIdx, meta, ...
        'Interactive', true);
    cfg.groups = groups_to_template(groups, meta.varNames);

    cfg.plotMode = prompt_plot_mode();
    if isempty(cfg.plotMode)
        cfg.cancelled = true;
        return;
    end

    cfg.saveIndividualScalograms = prompt_hsub_individual_scalograms();

    fprintf('Wavelet config locked: %d group(s), plot=%s, output → %s\n', ...
        numel(cfg.groups), cfg.plotMode, cfg.outRoot);
    for g = 1:numel(cfg.groups)
        fprintf('  Group "%s": %d mice\n', cfg.groups(g).name, numel(cfg.groups(g).colIdx));
    end

    up.lastWaveInput = cfg.rawPath;
    save_user_paths(up);
end
