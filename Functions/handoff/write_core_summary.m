function summaryPath = write_core_summary(handoffDir, fileStem, payload)
%WRITE_CORE_SUMMARY Script 2 → 3 handoff MAT (no full CWT arrays).

    ensure_dir(handoffDir);
    cfg = core_defaults();
    summaryPath = fullfile(handoffDir, [cfg.handoff.summaryPrefix fileStem '.mat']);

    payload.handoff_version = cfg.handoff.version;
    payload.fileStem = fileStem;
    payload.created = datetime('now');

    save(summaryPath, '-struct', 'payload', '-v7.3');
    fprintf('Handoff written: %s\n', summaryPath);
end
