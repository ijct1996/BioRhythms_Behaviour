function out = extended_script7_run(extHandoffDir, cfg)
%EXTENDED_SCRIPT7_RUN Phase events + profiles (local Extended implementation).

    if nargin < 2 || isempty(cfg)
        cfg = extended_defaults();
    end

    cfg = extended_apply_plot_cfg(cfg);
    out = extended_script7_profiles_run(extHandoffDir, cfg);
    fprintf('\nExtended Script 7 complete.\n');
end
