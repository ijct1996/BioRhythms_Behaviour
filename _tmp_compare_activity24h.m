legacyX = "c:/Users/User/OneDrive/Desktop/Files/Kent/05. MATLAB/03. Mice_behav/Codes/Final Versions/Pub_Res/SelectedValidatedUR_PublicationProfiles/Tables/SelectedValidatedUR_PublicationProfiles_Output.xlsx";
curX = "C:/Users/User/Dev/Cursor/Research/Chronobiology/Data/BioRhythms Behaviour Data/Results/C57_LP/ExtendedHandoff/SelectedValidatedUR_PublicationProfiles/Tables/SelectedValidatedUR_PublicationProfiles_Output.xlsx";
L = readtable(legacyX, "Sheet", "ActivityComponent_24h", "VariableNamingRule", "preserve");
C = readtable(curX, "Sheet", "ActivityComponent_24h", "VariableNamingRule", "preserve");
fprintf("LEGACY rows=%d\n", height(L));
fprintf("CURRENT rows=%d\n", height(C));
Ll = L(L.Photoperiod_h==12 & string(L.BandName)=="UR_1_3" & contains(string(L.ClusterID),"C01"), :);
Cc = C(C.Photoperiod_h==12 & string(C.BandName)=="UR_1_3" & contains(string(C.ClusterID),"C01"), :);
fprintf("LEGACY subset=%d CURRENT subset=%d\n", height(Ll), height(Cc));
GmL = groupsummary(Ll, "ZTBinCenter_h", "mean", "Activity_zscored");
GmC = groupsummary(Cc, "ZTBinCenter_h", "mean", "Activity_zscored");
fprintf("LEGACY min=%.3f max=%.3f\n", min(GmL.mean_Activity_zscored), max(GmL.mean_Activity_zscored));
fprintf("CURRENT min=%.3f max=%.3f\n", min(GmC.mean_Activity_zscored), max(GmC.mean_Activity_zscored));
i = @(G,z) mean(G.mean_Activity_zscored(abs(G.ZTBinCenter_h-z)<0.35), "omitnan");
fprintf("LEGACY @1=%.3f @12=%.3f @18=%.3f\n", i(GmL,1), i(GmL,12), i(GmL,18));
fprintf("CURRENT @1=%.3f @12=%.3f @18=%.3f\n", i(GmC,1), i(GmC,12), i(GmC,18));
fprintf("LEGACY frac0=%.3f CURRENT frac0=%.3f\n", mean(abs(GmL.mean_Activity_zscored)<0.05), mean(abs(GmC.mean_Activity_zscored)<0.05));
CLS = readtable(legacyX, "Sheet", "ClusterSummary", "VariableNamingRule", "preserve");
CSC = readtable(curX, "Sheet", "ClusterSummary", "VariableNamingRule", "preserve");
disp("LEGACY clusters"); disp(CLS);
disp("CURRENT clusters"); disp(CSC);
try
W=readtable(curX,"Sheet","ActivityWarnings","VariableNamingRule","preserve");
fprintf("Warnings=%d\n", height(W));
if height(W)>0, disp(W(1:min(12,height(W)),:)); end
catch, fprintf("No warnings\n"); end