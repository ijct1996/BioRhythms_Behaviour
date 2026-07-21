function out = extended_script5_run(extHandoffDir, cfg)
%EXTENDED_SCRIPT5_RUN Transition resync + FDR + photoperiod gradient (Script 5).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    fprintf('\n=== Extended Script 5: Transition resync + FDR ===\n');
    fprintf('Plot mode: %s | Primary bands: %s\n', cfg.plotMode, strjoin(cfg.bands.primaryUR, ', '));

    mapPath = '';
    if nargin >= 1 && ~isempty(extHandoffDir)
        cand = fullfile(extHandoffDir, 'RawVsSelectiveHSub_PeriodValidation', ...
            'HSubSupported_PeriodMap.mat');
        if isfile(cand)
            mapPath = cand;
            handoffDir = extHandoffDir;
        elseif isfile(fullfile(extHandoffDir, 'AcrossPhotoperiod_Input', ...
                'RawVsSelectiveHSub_PeriodValidation', 'HSubSupported_PeriodMap.mat'))
            handoffDir = fullfile(extHandoffDir, 'AcrossPhotoperiod_Input');
            mapPath = fullfile(handoffDir, 'RawVsSelectiveHSub_PeriodValidation', ...
                'HSubSupported_PeriodMap.mat');
        else
            handoffDir = extHandoffDir;
        end
    else
        handoffDir = '';
    end

    opts = struct('cfg', cfg);
    if ~isempty(mapPath)
        opts.mapPath = mapPath;
        opts.handoffDir = handoffDir;
        resyncAll = extended_ridge_resync_run(opts);
    else
        resyncAll = extended_ridge_resync_run(handoffDir, cfg);
    end

    out = struct();
    out.allDays = resyncAll;

    if ~extended_script5_resync_succeeded_(resyncAll)
        warning('extended_script5_run:MainPassFailed', ...
            'Main transition-resync pass did not complete; skipping gradient plots and stable-days sensitivity.');
        fprintf('\nExtended Script 5 aborted (main pass incomplete).\n');
        return;
    end

    extended_plot_transition_gradient(resyncAll, cfg);

    if cfg.plot.generateFigures && isfield(resyncAll, 'outXLSX')
        try
            extended_plot_ridge_power_takeaway(resyncAll.outXLSX, cfg);
        catch ME
            warning('extended_script5_run:TakeawayFailed', '%s', ME.message);
        end
    end

    % Reuse paths resolved by the main pass (avoids re-prompting on sensitivity pass)
    [mapPath, handoffDir] = extended_script5_resolved_paths_(mapPath, handoffDir, resyncAll);

    % Stable-days sensitivity (optional second pass)
    if cfg.ridge.runStableDaysSensitivity
        fprintf('Running stable-days sensitivity (Time_days >= %.1f)...\n', cfg.ridge.stableDaysMin);
        if isempty(mapPath) || isempty(handoffDir)
            warning('extended_script5_run:StableDaysSkipped', ...
                ['Stable-days sensitivity skipped: could not reuse map/handoff paths from main pass.\n', ...
                'Main-pass outputs remain valid at: %s'], resyncAll.outRoot);
        else
            cfgStable = cfg;
            cfgStable.ridge.excludeInitialDays = true;
            cfgStable.ridge.minDayForStableResync = cfg.ridge.stableDaysMin;
            optsStable = struct( ...
                'cfg', cfgStable, ...
                'mapPath', mapPath, ...
                'handoffDir', handoffDir, ...
                'outputSubdir', 'StableDaysSensitivity');
            try
                resyncStable = extended_ridge_resync_run(optsStable);
            catch ME
                warning('extended_script5_run:StableDaysFailed', ...
                    'Stable-days sensitivity pass failed: %s\nMain-pass outputs remain valid at: %s', ...
                    ME.message, resyncAll.outRoot);
                resyncStable = struct();
            end
            if extended_script5_resync_succeeded_(resyncStable)
                out.stableDays = resyncStable;
                extended_plot_transition_gradient(resyncStable, cfg);
                fprintf('  Stable-days output: %s\n', resyncStable.outRoot);
            else
                warning('extended_script5_run:StableDaysSkipped', ...
                    'Stable-days sensitivity pass did not complete.\nMain-pass outputs remain valid at: %s', ...
                    resyncAll.outRoot);
            end
        end
    end

    fprintf('\nExtended Script 5 complete.\n');
    if isfield(resyncAll, 'outRoot')
        fprintf('  Main output: %s\n', resyncAll.outRoot);
    end
end

function ok = extended_script5_resync_succeeded_(resync)
    ok = isstruct(resync) && isfield(resync, 'outRoot') && ~isempty(resync.outRoot) ...
        && isfield(resync, 'outXLSX') && isfile(resync.outXLSX);
end

function [mapPath, handoffDir] = extended_script5_resolved_paths_(mapPath, handoffDir, resyncAll)
    if isfield(resyncAll, 'mapPath') && ~isempty(resyncAll.mapPath)
        mapPath = char(string(resyncAll.mapPath));
    end
    if isfield(resyncAll, 'handoffDir') && ~isempty(resyncAll.handoffDir)
        handoffDir = char(string(resyncAll.handoffDir));
    elseif isfield(resyncAll, 'outRoot') && ~isempty(resyncAll.outRoot)
        handoffDir = fileparts(resyncAll.outRoot);
        if endsWith(handoffDir, 'StableDaysSensitivity')
            handoffDir = fileparts(handoffDir);
        end
        if endsWith(handoffDir, 'Ultradian_RidgePhase_Resync')
            handoffDir = fileparts(handoffDir);
        end
    end
    if isempty(mapPath) && ~isempty(handoffDir)
        cand = fullfile(handoffDir, 'RawVsSelectiveHSub_PeriodValidation', ...
            'HSubSupported_PeriodMap.mat');
        if isfile(cand)
            mapPath = cand;
        end
    end
end
