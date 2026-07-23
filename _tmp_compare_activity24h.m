legacyX = 'c:\Users\User\OneDrive\Desktop\Files\Kent\05. MATLAB\03. Mice_behav\Codes\Final Versions\Pub_Res\SelectedValidatedUR_PublicationProfiles\Tables\SelectedValidatedUR_PublicationProfiles_Output.xlsx';
curX = 'C:\Users\User\Dev\Cursor\Research\Chronobiology\Data\BioRhythms Behaviour Data\Results\C57_LP\ExtendedHandoff\SelectedValidatedUR_PublicationProfiles\Tables\SelectedValidatedUR_PublicationProfiles_Output.xlsx';

L = readtable(legacyX, 'Sheet', 'ActivityComponent_24h', 'VariableNamingRule', 'preserve');
C = readtable(curX, 'Sheet', 'ActivityComponent_24h', 'VariableNamingRule', 'preserve');

fprintf('LEGACY rows=%d\n', height(L));
fprintf('CURRENT rows=%d\n', height(C));
disp(L.Properties.VariableNames);
disp(C.Properties.VariableNames);

pick = @(T) T(T.Photoperiod_h == 12 & string(T.BandName) == "UR_1_3" & contains(string(T.ClusterID), "C01"), :);
Ll = pick(L);
Cc = pick(C);

if ismember('File', Ll.Properties.VariableNames)
    nMiceL = numel(unique(string(Ll.File) + "|" + string(Ll.SignalID)));
else
    nMiceL = numel(unique(string(Ll.SignalID)));
end
nMiceC = numel(unique(string(Cc.File) + "|" + string(Cc.SignalID)));
fprintf('LEGACY L12 UR13 C01 rows=%d mice=%d\n', height(Ll), nMiceL);
fprintf('CURRENT L12 UR13 C01 rows=%d mice=%d\n', height(Cc), nMiceC);

GmL = groupsummary(Ll, 'ZTBinCenter_h', 'mean', 'Activity_zscored');
GmC = groupsummary(Cc, 'ZTBinCenter_h', 'mean', 'Activity_zscored');
fprintf('LEGACY ZT bins=%d mean act min=%.3f max=%.3f\n', height(GmL), min(GmL.mean_Activity_zscored), max(GmL.mean_Activity_zscored));
fprintf('CURRENT ZT bins=%d mean act min=%.3f max=%.3f\n', height(GmC), min(GmC.mean_Activity_zscored), max(GmC.mean_Activity_zscored));

near = @(G, z) mean(G.mean_Activity_zscored(abs(G.ZTBinCenter_h - z) < 0.35), 'omitnan');
fprintf('LEGACY mean@~1h=%.3f @~12h=%.3f @~18h=%.3f\n', near(GmL, 1), near(GmL, 12), near(GmL, 18));
fprintf('CURRENT mean@~1h=%.3f @~12h=%.3f @~18h=%.3f\n', near(GmC, 1), near(GmC, 12), near(GmC, 18));

% per-bin SD across mice at ZT~1 and ZT~12 to see if individuals collapsed
sd_at = @(T, z) std(T.Activity_zscored(abs(T.ZTBinCenter_h - z) < 0.35), 'omitnan');
fprintf('LEGACY SD@~1h=%.3f @~12h=%.3f\n', sd_at(Ll, 1), sd_at(Ll, 12));
fprintf('CURRENT SD@~1h=%.3f @~12h=%.3f\n', sd_at(Cc, 1), sd_at(Cc, 12));

CLS = readtable(legacyX, 'Sheet', 'ClusterSummary', 'VariableNamingRule', 'preserve');
CSC = readtable(curX, 'Sheet', 'ClusterSummary', 'VariableNamingRule', 'preserve');
keep = {'ClusterID', 'BandName', 'ClusterRank', 'PeriodCentre_h', 'PeriodLow_h', 'PeriodHigh_h', 'CandidateCount'};
disp('LEGACY ClusterSummary:');
disp(CLS(:, intersect(CLS.Properties.VariableNames, keep, 'stable')));
disp('CURRENT ClusterSummary:');
disp(CSC(:, intersect(CSC.Properties.VariableNames, keep, 'stable')));

try
    W = readtable(curX, 'Sheet', 'ActivityWarnings', 'VariableNamingRule', 'preserve');
    fprintf('CURRENT ActivityWarnings rows=%d\n', height(W));
    if height(W) > 0
        disp(W(1:min(15, height(W)), :));
    end
catch
    fprintf('No ActivityWarnings\n');
end

% Filter window columns if present
if all(ismember({'FilterLow_h', 'FilterHigh_h'}, Cc.Properties.VariableNames))
    uC = unique(Cc(:, {'FilterLow_h', 'FilterHigh_h'}));
    disp('CURRENT filter windows:');
    disp(uC);
end
if all(ismember({'FilterLow_h', 'FilterHigh_h'}, Ll.Properties.VariableNames))
    uL = unique(Ll(:, {'FilterLow_h', 'FilterHigh_h'}));
    disp('LEGACY filter windows:');
    disp(uL);
end
