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
    extended_plot_transition_gradient(resyncAll, cfg);

    if cfg.plot.generateFigures && isfield(resyncAll, 'outXLSX')
        try
            extended_plot_ridge_power_takeaway(resyncAll.outXLSX, cfg);
        catch ME
            warning('extended_script5_run:TakeawayFailed', '%s', ME.message);
        end
    end

    % Stable-days sensitivity (optional second pass)
    if cfg.ridge.runStableDaysSensitivity
        fprintf('Running stable-days sensitivity (Time_days >= %.1f)...\n', cfg.ridge.stableDaysMin);
        cfgStable = cfg;
        cfgStable.ridge.excludeInitialDays = true;
        cfgStable.ridge.minDayForStableResync = cfg.ridge.stableDaysMin;
        if ~isempty(mapPath)
            optsStable = struct('cfg', cfgStable, 'mapPath', mapPath, 'handoffDir', handoffDir);
            resyncStable = extended_ridge_resync_run(optsStable);
        else
            resyncStable = extended_ridge_resync_run(handoffDir, cfgStable);
        end
        out.stableDays = resyncStable;
        if isfield(resyncStable, 'outRoot')
            stableDir = fullfile(resyncStable.outRoot, 'StableDaysSensitivity');
            extended_period_gate_ensure_dir(stableDir);
            extended_plot_transition_gradient(resyncStable, cfg);
        end
    end

    fprintf('\nExtended Script 5 complete.\n');
    if isfield(resyncAll, 'outRoot')
        fprintf('  Output: %s\n', resyncAll.outRoot);
    end
end
