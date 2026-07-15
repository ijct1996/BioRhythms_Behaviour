function pal = collaborator_palette()
%COLLABORATOR_PALETTE Colours from Figure Graph Formatting Information_V2.pptx
%
%   pal = collaborator_palette()

    pal.C57_grouped = [64, 128, 128] / 255;       % #408080, bars 50% transparency applied at plot time
    pal.C57_alpha = 0.5;

    pal.NR2B_wt     = [0, 0, 4] / 255;            % #000004 +/+
    pal.NR2B_ko     = [59, 15, 112] / 255;        % #3B0F70 −/−

    pal.NR2B_male   = [200, 222, 229] / 255;      % #C8DEE5 — KO / lighter male
    pal.NR2B_female = [238, 204, 249] / 255;      % #EECCF9 — KO / lighter female

    % Sex × genotype (4-group plots): Ctrl = darker shade, KO = lighter (slide colours above)
    pal.NR2B_male_ctrl   = blend_rgb(pal.NR2B_wt, pal.NR2B_male, 0.55);   % dark blue-grey
    pal.NR2B_female_ctrl = blend_rgb(pal.NR2B_ko, pal.NR2B_female, 0.45); % dark pink-purple

    pal.phase_a     = [173, 7, 227] / 255;        % #AD07E3
    pal.phase_b     = [9, 153, 99] / 255;         % #099963

    pal.axes.fontName = 'Calibri';
    pal.axes.fontSize = 12;
    pal.axes.fontWeight = 'bold';
    pal.axes.lineWidth = 1;

    pal.publication.fontName = 'Times New Roman';
end

function c = blend_rgb(dark, light, darkWeight)
    w = min(max(darkWeight, 0), 1);
    c = w * dark(:)' + (1 - w) * light(:)';
end
