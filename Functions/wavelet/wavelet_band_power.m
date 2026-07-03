function bandPower = wavelet_band_power(powerSpec, periods_hours, cfg)
%WAVELET_BAND_POWER Integrate mean power in CR/UR bands for co-expression.

    if nargin < 3, cfg = core_defaults(); end
    bandNames = cfg.bands.coexpression;
    bandPower = struct();

    for i = 1:numel(bandNames)
        bn = bandNames{i};
        lim = cfg.bands.(bn);
        mask = periods_hours >= lim(1) & periods_hours <= lim(2);
        if any(mask)
            bandPower.(bn) = mean(powerSpec(mask), 'omitnan');
        else
            bandPower.(bn) = NaN;
        end
    end
end
