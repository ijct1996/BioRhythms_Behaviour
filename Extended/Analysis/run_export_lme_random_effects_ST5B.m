function run_export_lme_random_effects_ST5B(lmeOutputsXlsx, outDir)
%RUN_EXPORT_LME_RANDOM_EFFECTS_ST5B Independent ST5B random-effects export.
%
%   Refits the four AdditiveOnly manuscript LMEs and writes mouse / file /
%   residual variance components for Supplementary Table 5B.
%
%   Does NOT re-run Script 6. Requires Script 6 tables already present:
%     AcrossPhotoperiod_LME/Tables/AcrossPhotoperiod_Outputs.xlsx
%
%   Usage:
%     run_export_lme_random_effects_ST5B
%     run_export_lme_random_effects_ST5B(xlsxPath, outDir)

    setup_extended_paths();

    if nargin < 1
        lmeOutputsXlsx = [];
    end
    if nargin < 2
        outDir = [];
    end

    out = extended_export_lme_random_effects_ST5B(lmeOutputsXlsx, outDir);
    fprintf('\nST5B export complete (%d rows).\n', height(out.rows));
    disp(out.rows(:, {'Outcome', 'RandomEffect', 'VarianceComponent'}));
end
