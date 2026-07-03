function ratioTable = across_compute_coexpression(bandTable, cfg)
%ACROSS_COMPUTE_COEXPRESSION CR vs UR band ratios per mouse.

    if nargin < 2, cfg = core_defaults(); end

    urBands = {'UR_1_3', 'UR_3_6', 'UR_6_12', 'UR_12_18'};
    ratioTable = bandTable(:, {'SignalID', 'Group', 'LightDuration_h', 'CR_20_28'});

    cr = bandTable.CR_20_28;
    for i = 1:numel(urBands)
        ur = bandTable.(urBands{i});
        ratioTable.(['CR_to_' urBands{i}]) = cr ./ (ur + eps);
        ratioTable.([urBands{i} '_fraction_of_UR_total']) = ur ./ ...
            (bandTable.UR_1_3 + bandTable.UR_3_6 + bandTable.UR_6_12 + bandTable.UR_12_18 + eps);
    end

    ratioTable.UR_total = bandTable.UR_1_3 + bandTable.UR_3_6 + bandTable.UR_6_12 + bandTable.UR_12_18;
    ratioTable.CR_to_UR_total = cr ./ (ratioTable.UR_total + eps);
end
