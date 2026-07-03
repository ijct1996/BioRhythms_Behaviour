%% RUN_ACROSS_PHOTOPERIOD — Core v1 Script 3
% Reads Handoff/ + Script 1 residuals. RAW (exploratory) and HSub-validated UR
% (confirmatory) co-expression / period comparison; stitched scalograms.
%
%   setup_paths();
%   run_across_photoperiod

function run_across_photoperiod()
    setup_paths();
    cfg = core_defaults();

    fprintf('BioRhythms Core v1 — Across photoperiod (%s)\n', cfg.version);

    plotMode = prompt_plot_mode();
    if isempty(plotMode)
        fprintf('Cancelled.\n');
        return;
    end

    handoffDir = uigetdir(pwd, 'Select cohort Handoff/ folder');
    if isequal(handoffDir, 0), return; end

    listStr = cfg.coreCohorts;
    [sel, ok] = listdlg('PromptString', 'Select cohort tag:', ...
        'SelectionMode', 'single', 'ListString', listStr);
    if ~ok, return; end
    cohortTag = listStr{sel};

    outRoot = uigetdir(handoffDir, 'Select cohort results root (parent of 03_AcrossPhotoperiod)');
    if isequal(outRoot, 0), return; end

    across_run_cohort(handoffDir, outRoot, cohortTag, plotMode);
    fprintf('Done.\n');
end
