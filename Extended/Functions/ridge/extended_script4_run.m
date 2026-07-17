function out = extended_script4_run(coreHandoffDir, cfg)
%EXTENDED_SCRIPT4_RUN Ridge handoff from Core + CarryForward gate (Extended Script 4).
%
%   out = extended_script4_run()
%   out = extended_script4_run(coreHandoffDir)
%   out = extended_script4_run(coreHandoffDir, cfg)
%
%   If WP_Summary__*.mat already exist under ExtendedHandoff, offers
%   "Gate only" so a successful ridge handoff need not be recomputed.

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    fprintf('\n=== Extended Script 4: Ridge handoff + CarryForward validation ===\n');
    fprintf('Plot mode: %s\n', cfg.plotMode);

    if nargin < 1 || isempty(coreHandoffDir)
        coreHandoffDir = uigetdir(pwd, 'Select Core cohort Handoff/ folder (CoreSummary__*.mat)');
        if isequal(coreHandoffDir, 0)
            fprintf('Cancelled.\n');
            out = struct();
            return;
        end
    end

    entries = across_load_handoff_cohort(coreHandoffDir);
    extRoot = fullfile(fileparts(coreHandoffDir), 'ExtendedHandoff');
    extHandoffDir = fullfile(extRoot, 'AcrossPhotoperiod_Input');
    extended_period_gate_ensure_dir(extHandoffDir);

    existing = dir(fullfile(extHandoffDir, 'WP_Summary__*.mat'));
    recompute = true;
    if ~isempty(existing)
        choice = questdlg(sprintf( ...
            ['Found %d existing WP_Summary__*.mat under ExtendedHandoff.\n\n' ...
             'Gate only  — run CarryForward on existing ridge handoff (recommended if ridge already finished).\n' ...
             'Recompute  — rebuild ridge handoff for all photoperiods (slow).\n'], ...
            numel(existing)), ...
            'Script 4 resume', ...
            'Gate only', 'Recompute all', 'Cancel', 'Gate only');
        if isempty(choice) || strcmp(choice, 'Cancel')
            fprintf('Cancelled.\n');
            out = struct();
            return;
        end
        recompute = strcmp(choice, 'Recompute all');
    end

    if recompute
        ridgeSummary = extended_ridge_handoff_from_core(entries, extHandoffDir, cfg);
    else
        fprintf('Reusing existing ridge handoff (%d WP_Summary files).\n', numel(existing));
        ridgeSummary = struct();
        ridgeSummary.extHandoffDir = extHandoffDir;
        ridgeSummary.reused = true;
        ridgeSummary.nExisting = numel(existing);
    end

    validationMap = extended_period_gate_run(extHandoffDir, cfg);

    out = struct();
    out.coreHandoffDir = coreHandoffDir;
    out.extHandoffDir = extHandoffDir;
    out.ridgeSummary = ridgeSummary;
    out.validationMap = validationMap;
    out.mapPath = fullfile(extHandoffDir, 'RawVsSelectiveHSub_PeriodValidation', ...
        'HSubSupported_PeriodMap.mat');

    fprintf('\nExtended Script 4 complete.\n');
    fprintf('  Extended handoff: %s\n', extHandoffDir);
    fprintf('  CarryForward map: %s\n', out.mapPath);
end
