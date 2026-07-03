function cfg = prompt_hsub_batch_config()
%PROMPT_HSUB_BATCH_CONFIG Collect all Script 1 settings upfront (unsupervised batch).

    up = load_user_paths();
    startDir = up.lastHSubInput;
    if isempty(startDir) || ~isfolder(startDir)
        startDir = pwd;
    end

    cfg = struct();
    cfg.version = core_defaults().version;

    [cfg.files, cfg.inPath, cfg.inputMode] = prompt_input_files(startDir, ...
        'Select Pipeline_Input xlsx file(s)');
    if isempty(cfg.files)
        cfg.cancelled = true;
        return;
    end
    cfg.cancelled = false;
    fprintf('HSub batch: %s mode, %d file(s)\n', cfg.inputMode, numel(cfg.files));

    cfg.outRoot = uigetdir(cfg.inPath, ...
        'Select/create 01_HarmonicSubtraction output folder');
    if isequal(cfg.outRoot, 0)
        cfg.cancelled = true;
        return;
    end
    ensure_dir(cfg.outRoot);

    firstFile = fullfile(cfg.inPath, cfg.files{1});
    [~, meta] = read_behav_excel(firstFile);
    cfg.mouseNames = prompt_mouse_columns_shared(meta.varNames, meta.mouseIdx);
    if isempty(cfg.mouseNames)
        cfg.cancelled = true;
        return;
    end

    cfg.groups = prompt_hsub_groups(meta.varNames, meta.mouseIdx, meta);
    if isempty(cfg.groups)
        cfg.cancelled = true;
        return;
    end

    cfg.plotMode = prompt_plot_mode();
    if isempty(cfg.plotMode)
        cfg.cancelled = true;
        return;
    end

    cfg.saveIndividualScalograms = prompt_hsub_individual_scalograms();

    fprintf('HSub config locked: %d mice, %d group(s), plot=%s, output → %s\n', ...
        numel(cfg.mouseNames), numel(cfg.groups), cfg.plotMode, cfg.outRoot);

    up.lastHSubInput = cfg.inPath;
    save_user_paths(up);
end
