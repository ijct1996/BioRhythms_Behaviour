function peaks = wavelet_find_period_peaks(powerSpec, periods_hours, topN)
%WAVELET_FIND_PERIOD_PEAKS Top-N peaks from mean log10 power spectrum.

    if nargin < 3, topN = core_defaults().wavelet.topNPeaks; end
    logP = log10(powerSpec + eps);
    [pks, locs, w, prom] = findpeaks(logP, 'MinPeakProminence', 0);

    if isempty(pks)
        peaks = table();
        return;
    end

    [~, ord] = sort(pks, 'descend');
    ord = ord(1:min(topN, numel(ord)));
    peaks = table((1:numel(ord))', periods_hours(locs(ord)), pks(ord), prom(ord), w(ord), ...
        'VariableNames', {'PeakRank', 'PeakPeriod_hr', 'PeakValue_log10', 'PeakProminence', 'PeakWidth'});
end
