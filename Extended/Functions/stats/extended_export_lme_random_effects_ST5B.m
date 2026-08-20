function out = extended_export_lme_random_effects_ST5B(lmeOutputsXlsx, outDir)
%EXTENDED_EXPORT_LME_RANDOM_EFFECTS_ST5B Export AdditiveOnly LME VarCorr for ST5B.
%
%   Independent of Script 6 re-run: reads Script 6 tables and refits the four
%   manuscript AdditiveOnly models (Delta + Power × UR_1_3 + UR_3_6), then
%   writes mouse / file / residual variance components for Supp Table 5B.
%
%   out = extended_export_lme_random_effects_ST5B()
%   out = extended_export_lme_random_effects_ST5B(lmeOutputsXlsx, outDir)
%
%   Default input:
%     .../ExtendedHandoff/AcrossPhotoperiod_Input/AcrossPhotoperiod_LME/Tables/AcrossPhotoperiod_Outputs.xlsx
%   Default output:
%     .../Results/C57_LP/Supplemental_Tables_Publication/

    if nargin < 1 || isempty(lmeOutputsXlsx)
        lmeOutputsXlsx = default_outputs_xlsx_();
    end
    if nargin < 2 || isempty(outDir)
        outDir = default_outdir_();
    end
    lmeOutputsXlsx = char(string(lmeOutputsXlsx));
    outDir = char(string(outDir));
    if ~isfile(lmeOutputsXlsx)
        error('extended_export_lme_random_effects_ST5B:MissingXlsx', ...
            'AcrossPhotoperiod_Outputs.xlsx not found:\n%s\nRun Extended Script 6 once first.', lmeOutputsXlsx);
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fprintf('\n=== Export LME random effects (ST5B) ===\n');
    fprintf('Input:  %s\n', lmeOutputsXlsx);
    fprintf('Output: %s\n', outDir);

    Tpair = readtable(lmeOutputsXlsx, 'Sheet', 'CR_UR_Pairs_PerMouse');
    TBCS = readtable(lmeOutputsXlsx, 'Sheet', 'BandCondSummary_AnalysisUsed');

    urBands = ["UR_1_3"; "UR_3_6"];
    outcomeLabel = dictionary( ...
        "Delta|UR_1_3", "Δ UR 1–3 h", ...
        "Delta|UR_3_6", "Δ UR 3–6 h", ...
        "Power|UR_1_3", "Power UR 1–3 h", ...
        "Power|UR_3_6", "Power UR 3–6 h");

    rows = table();
    models = {};

    % --- Delta models (CR–UR pairs) ---
    Tpair = prepare_pair_table_(Tpair);
    for i = 1:numel(urBands)
        ub = urBands(i);
        D = Tpair(Tpair.UR_Band == categorical(ub), :);
        D = drop_unknown_sex_(D);
        fml = 'Delta_log10 ~ Photoperiod_h + Phase + Sex + (1|SignalID) + (1|File)';
        fprintf('Fitting Delta %s (n=%d) ...\n', ub, height(D));
        mdl = fitlme(D, fml);
        models{end+1} = struct('Metric', "Delta", 'Band', ub, 'Model', mdl, 'Formula', string(fml)); %#ok<AGROW>
        rows = [rows; variance_rows_(mdl, outcomeLabel("Delta|" + ub), "Delta", ub, fml)]; %#ok<AGROW>
    end

    % --- Power models (band power) ---
    TBP = prepare_power_table_(TBCS);
    for i = 1:numel(urBands)
        ub = urBands(i);
        D = TBP(TBP.BandName == categorical(ub), :);
        D = drop_unknown_sex_(D);
        fml = 'MeanBandPower_log10 ~ Photoperiod_h + Phase + Sex + (1|SignalID) + (1|File)';
        fprintf('Fitting Power %s (n=%d) ...\n', ub, height(D));
        mdl = fitlme(D, fml);
        models{end+1} = struct('Metric', "Power", 'Band', ub, 'Model', mdl, 'Formula', string(fml)); %#ok<AGROW>
        rows = [rows; variance_rows_(mdl, outcomeLabel("Power|" + ub), "Power", ub, fml)]; %#ok<AGROW>
    end

    % Manuscript display labels: SignalID ≡ mouse (matches Fig 3B legend).
    rows.RandomEffect = rename_re_for_manuscript_(rows.RandomEffect);

    csvPath = fullfile(outDir, 'SuppTable05B_LME_RandomEffects.csv');
    writetable(rows(:, {'Outcome', 'RandomEffect', 'VarianceComponent'}), csvPath);
    % Keep audit columns in a second file (GroupingVariable retains raw SignalID)
    auditPath = fullfile(outDir, 'SuppTable05B_LME_RandomEffects_Audit.csv');
    writetable(rows, auditPath);

    xlsxPath = fullfile(outDir, 'SuppTable05B_LME_RandomEffects.xlsx');
    if isfile(xlsxPath), delete(xlsxPath); end
    writetable(rows(:, {'Outcome', 'RandomEffect', 'VarianceComponent'}), xlsxPath, 'Sheet', 'ST5B');
    writetable(rows, xlsxPath, 'Sheet', 'Audit');

    footnotePath = fullfile(outDir, 'SuppTable05B_Footnote.txt');
    fid = fopen(footnotePath, 'w');
    fprintf(fid, ['Footnote: Random-effect labels match the Figure 3B / Methods formula ', ...
        'response ~ Photoperiod_h + Phase + Sex + (1|mouse) + (1|file). ', ...
        'In the fitted model the mouse grouping variable is coded as SignalID ', ...
        '(one SignalID = one mouse); File is the recording/photoperiod-block identifier. ', ...
        'Values are variance components from covarianceParameters(fitlme).\n']);
    fclose(fid);

    fprintf('Wrote:\n  %s\n  %s\n  %s\n  %s\n', csvPath, auditPath, xlsxPath, footnotePath);

    out = struct();
    out.rows = rows;
    out.csvPath = string(csvPath);
    out.auditPath = string(auditPath);
    out.xlsxPath = string(xlsxPath);
    out.models = models;
end

%% ---- helpers -----------------------------------------------------------

function T = prepare_pair_table_(T)
    T.Photoperiod_h = double(T.Photoperiod_h);
    T.SignalID = categorical(string(T.SignalID));
    T.File = categorical(string(T.File));
    T.Phase = categorical(string(T.Phase));
    T.UR_Band = categorical(string(T.UR_Band));
    if ~ismember('Sex', T.Properties.VariableNames)
        T = add_sex_column_(T);
    else
        T.Sex = categorical(string(T.Sex));
    end
end

function T = prepare_power_table_(BCS)
    BCS = BCS(string(BCS.Source) == "Raw", :);
    T = BCS(:, {'File', 'SignalID', 'Photoperiod_h', 'Phase', 'BandName', 'MeanBandPower_log10'});
    T.Photoperiod_h = double(T.Photoperiod_h);
    T.SignalID = categorical(string(T.SignalID));
    T.File = categorical(string(T.File));
    T.Phase = categorical(string(T.Phase));
    T.BandName = categorical(string(T.BandName));
    T = add_sex_column_(T);
end

function T = add_sex_column_(T)
    sid = string(T.SignalID);
    sx = strings(size(sid));
    for i = 1:numel(sid)
        s = upper(char(sid(i)));
        if contains(s, '_F') || startsWith(s, 'F_') || contains(s, '-F')
            sx(i) = "Female";
        elseif contains(s, '_M') || startsWith(s, 'M_') || contains(s, '-M')
            sx(i) = "Male";
        else
            sx(i) = "Unknown";
        end
    end
    T.Sex = categorical(sx);
end

function D = drop_unknown_sex_(D)
    sx = string(D.Sex);
    D = D(~ismissing(sx) & sx ~= "Unknown", :);
end

function rows = variance_rows_(mdl, outcome, metric, band, fml)
    [psi, mse] = covarianceParameters(mdl);
    gNames = grouping_names_(mdl, numel(psi));

    Outcome = strings(0, 1);
    RandomEffect = strings(0, 1);
    VarianceComponent = zeros(0, 1);
    MetricClass = strings(0, 1);
    BandName = strings(0, 1);
    Formula = strings(0, 1);
    GroupingVariable = strings(0, 1);

    for k = 1:numel(psi)
        v = psi{k};
        if isempty(v)
            continue;
        end
        % intercept-only RE → scalar variance
        if isscalar(v)
            varHat = double(v);
        else
            varHat = double(v(1, 1));
        end
        Outcome(end+1, 1) = string(outcome); %#ok<AGROW>
        RandomEffect(end+1, 1) = string(gNames{k}); %#ok<AGROW>
        VarianceComponent(end+1, 1) = varHat; %#ok<AGROW>
        MetricClass(end+1, 1) = string(metric); %#ok<AGROW>
        BandName(end+1, 1) = string(band); %#ok<AGROW>
        Formula(end+1, 1) = string(fml); %#ok<AGROW>
        GroupingVariable(end+1, 1) = string(gNames{k}); %#ok<AGROW>
    end

    Outcome(end+1, 1) = string(outcome);
    RandomEffect(end+1, 1) = "Residual";
    VarianceComponent(end+1, 1) = double(mse);
    MetricClass(end+1, 1) = string(metric);
    BandName(end+1, 1) = string(band);
    Formula(end+1, 1) = string(fml);
    GroupingVariable(end+1, 1) = "Residual";

    rows = table(Outcome, RandomEffect, VarianceComponent, MetricClass, BandName, GroupingVariable, Formula);
end

function names = grouping_names_(mdl, nPsi)
    % Prefer formula grouping names; fall back to Script 6 order.
    % Display rename to (1 | mouse) happens later in rename_re_for_manuscript_.
    names = cell(nPsi, 1);
    try
        g = mdl.Formula.GroupingVariableNames;
        if iscell(g) && numel(g) >= nPsi
            for i = 1:nPsi
                gi = g{i};
                if iscell(gi), gi = gi{1}; end
                names{i} = sprintf('(1 | %s)', char(string(gi)));
            end
            return;
        end
    catch
    end
    defaults = {'(1 | SignalID)', '(1 | File)'};
    for i = 1:nPsi
        if i <= numel(defaults)
            names{i} = defaults{i};
        else
            names{i} = sprintf('(1 | RE%d)', i);
        end
    end
end

function re = rename_re_for_manuscript_(re)
    re = string(re);
    re(re == "(1 | SignalID)" | re == "(1|SignalID)") = "(1 | mouse)";
    re(re == "(1 | File)" | re == "(1|File)") = "(1 | file)";
end

function p = default_outputs_xlsx_()
    p = fullfile( ...
        'C:', 'Users', 'User', 'Dev', 'Cursor', 'Research', 'Chronobiology', ...
        'Data', 'BioRhythms Behaviour Data', 'Results', 'C57_LP', ...
        'ExtendedHandoff', 'AcrossPhotoperiod_Input', 'AcrossPhotoperiod_LME', ...
        'Tables', 'AcrossPhotoperiod_Outputs.xlsx');
end

function p = default_outdir_()
    p = fullfile( ...
        'C:', 'Users', 'User', 'Dev', 'Cursor', 'Research', 'Chronobiology', ...
        'Data', 'BioRhythms Behaviour Data', 'Results', 'C57_LP', ...
        'Supplemental_Tables_Publication');
end
