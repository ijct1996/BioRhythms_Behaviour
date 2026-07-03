function style = pick_collaborator_group_style(grp, pal)
%PICK_COLLABORATOR_GROUP_STYLE Sex/genotype colours (collaborator slide spec).
%
%   Male   : #C8DEE5 (RGB 200,222,229), open circle, 1 pt black edge
%   Female : #EECCF9 (RGB 238,204,249), diamond, 1 pt black edge
%
%   Recognises group labels: F, M, Female, Male (case-insensitive).
%   Test female before male — contains('female','male') is true in MATLAB.

    if nargin < 2, pal = collaborator_palette(); end

    g = lower(strtrim(char(grp)));
    style = struct();
    style.lineColor = pal.NR2B_ko;
    style.marker = 'o';
    style.markerFaceColor = 'none';
    style.markerEdgeColor = [0 0 0];
    style.markerSize = 6;
    style.lineWidth = 1.5;

    if is_collaborator_female_group(g)
        style.lineColor = pal.NR2B_female;
        style.marker = 'd';
        style.markerFaceColor = pal.NR2B_female;
        style.markerEdgeColor = [0 0 0];
    elseif is_collaborator_male_group(g)
        style.lineColor = pal.NR2B_male;
        style.marker = 'o';
        style.markerFaceColor = 'none';
        style.markerEdgeColor = [0 0 0];
    elseif contains(g, 'c57') || strcmp(g, 'all')
        style.lineColor = pal.C57_grouped;
        style.marker = 'o';
        style.markerFaceColor = 'none';
        style.markerEdgeColor = [0 0 0];
    end
end

function tf = is_collaborator_female_group(g)
    tf = any(strcmp(g, {'f', 'female', 'fem', 'females'})) || contains(g, 'female');
end

function tf = is_collaborator_male_group(g)
    tf = any(strcmp(g, {'m', 'male', 'males'})) ...
        || (contains(g, 'male') && ~contains(g, 'female'));
end
