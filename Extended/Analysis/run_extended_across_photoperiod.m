function run_extended_across_photoperiod(varargin)
%RUN_EXTENDED_ACROSS_PHOTOPERIOD Deprecated — use run_extended_script6_across_lme.
    warning('run_extended_across_photoperiod:Deprecated', ...
        'Use run_extended_script6_across_lme instead.');
    setup_extended_paths();
    if nargin > 0
        run_extended_script6_across_lme(varargin{1});
    else
        run_extended_script6_across_lme();
    end
end
