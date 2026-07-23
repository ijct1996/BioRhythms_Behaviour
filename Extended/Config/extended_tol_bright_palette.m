function pal = extended_tol_bright_palette()
%EXTENDED_TOL_BRIGHT_PALETTE Paul Tol Bright palette + Script 8 semantic maps.
%
%   Colour contract: one encoded dimension per panel; roles do not overlap.
%   See Extended Script 8 Manifest for figure-specific usage.

    pal.name = 'PaulTolBright';
    pal.base = [
        68 119 170;
        102 204 238;
        34 136 51;
        204 187 68;
        238 102 119;
        170 51 119;
        187 187 187] / 255;

    pal.hex = {'#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB'};
    pal.labels = {'Blue', 'Cyan', 'Green', 'Yellow', 'Red', 'Purple', 'Grey'};

    %% Reserved roles (global)
    pal.dl = pal.base(1, :);       % DL lights-on
    pal.ld = pal.base(5, :);       % LD lights-off
    pal.cr = pal.base(7, :);       % CR band / dark-phase schematic
    pal.female = pal.base(3, :);   % Sex inset only
    pal.male = pal.base(4, :);     % Sex inset only
    pal.pooled = pal.base(7, :) * 0.55 + [1 1 1] * 0.45;

    pal.l12 = pal.base(1, :);      % 24h profile comparison
    pal.l24 = pal.base(6, :);      % LL endpoint (label/facet; line in 24h plots)

    %% Primary UR hero bands
    pal.band = containers.Map();
    pal.band('UR_1_3') = pal.base(2, :);
    pal.band('UR_3_6') = pal.base(6, :);
    pal.band('UR_6_9') = pal.base(1, :);
    pal.band('UR_9_12') = pal.base(5, :);
    pal.band('UR_12_18') = pal.base(3, :);
    pal.band('CR_20_28') = pal.base(7, :);

    pal.primaryUR = ["UR_1_3", "UR_3_6"];
    pal.allUR = ["UR_1_3", "UR_3_6", "UR_6_9", "UR_9_12", "UR_12_18"];

    pal.coherenceFacets = [12, 14, 16, 18, 20, 22];  % entrained photoperiods only (no LL)
    pal.coherenceYMax = 0.65;
    pal.coherenceXlim = [-6, 6];

    pal.diverging = [pal.base(1, :); 1 1 1; pal.base(5, :)];

    pal.fontName = 'Times New Roman';
    pal.titleWeight = 'bold';
    pal.labelWeight = 'bold';
    pal.axesLineWidth = 1.0;
    pal.tickDir = 'out';
end
