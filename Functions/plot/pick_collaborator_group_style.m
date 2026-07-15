function style = pick_collaborator_group_style(grp, pal)
%PICK_COLLABORATOR_GROUP_STYLE Sex / sex×genotype colours (collaborator spec).
%
%   2 groups (F, M):     female #EECCF9 diamond; male #C8DEE5 open circle
%   4 groups (F/M × Ctrl/KO): same sex hue; Ctrl darker, KO lighter
%
%   Labels: "F - Ctrl", "F - KO", "M - Ctrl", "M - KO" (and variants)

    if nargin < 2, pal = collaborator_palette(); end

    meta = parse_collaborator_group(grp);
    style = struct();
    style.lineColor = pal.NR2B_ko;
    style.marker = 'o';
    style.markerFaceColor = 'none';
    style.markerEdgeColor = [0 0 0];
    style.markerSize = 6;
    style.lineWidth = 1.5;

    if strcmp(meta.sex, 'F')
        style.marker = 'd';
        if strcmp(meta.genotype, 'Ctrl')
            style.lineColor = pal.NR2B_female_ctrl;
            style.markerFaceColor = pal.NR2B_female_ctrl;
        else
            style.lineColor = pal.NR2B_female;
            style.markerFaceColor = pal.NR2B_female;
        end
        style.markerEdgeColor = [0 0 0];
    elseif strcmp(meta.sex, 'M')
        style.marker = 'o';
        if strcmp(meta.genotype, 'Ctrl')
            style.lineColor = pal.NR2B_male_ctrl;
            style.markerFaceColor = 'none';
        else
            style.lineColor = pal.NR2B_male;
            style.markerFaceColor = 'none';
        end
        style.markerEdgeColor = [0 0 0];
    elseif contains(meta.raw, 'c57') || strcmp(meta.raw, 'all')
        style.lineColor = pal.C57_grouped;
        style.marker = 'o';
        style.markerFaceColor = 'none';
        style.markerEdgeColor = [0 0 0];
    end
end
