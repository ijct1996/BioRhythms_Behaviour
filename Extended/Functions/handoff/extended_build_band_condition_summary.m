function BCS = extended_build_band_condition_summary(periodCandidates)
%EXTENDED_BUILD_BAND_CONDITION_SUMMARY BandConditionSummary from Script 4 candidates.
%
%   Script 4 ridge handoff writes PeriodCandidates_Long but not the full Kent B
%   BandConditionSummary table. Script 6 builds the latter from Raw candidates.

    if isempty(periodCandidates) || height(periodCandidates) == 0
        BCS = table();
        return;
    end

    T = periodCandidates;
    T.Source = string(T.Source);
    T = T(strcmpi(T.Source, 'Raw'), :);
    if isempty(T)
        BCS = table();
        return;
    end

    n = height(T);
    BCS = table();
    BCS.File = string(T.File);
    BCS.SignalID = string(T.SignalID);
    BCS.Source = string(T.Source);
    BCS.Photoperiod_h = double(T.Photoperiod_h);
    BCS.LightStateValue = repmat("NA", n, 1);
    BCS.Phase = string(T.Phase);
    BCS.BandName = string(T.BandName);
    BCS.MeanBandPower_log10 = double(T.MeanBandPower_log10);
    BCS.SDBandPower_log10 = double(T.SDBandPower_log10);
    BCS.MeanBandPower_linear = 10 .^ BCS.MeanBandPower_log10;
    BCS.SDBandPower_linear = nan(n, 1);
    BCS.FracBandPower_linear = nan(n, 1);
    if ismember('MedianRidgePeriod_h', T.Properties.VariableNames)
        BCS.MeanRidgePeriod = double(T.MedianRidgePeriod_h);
    else
        BCS.MeanRidgePeriod = nan(n, 1);
    end
    if ismember('IQR_RidgePeriod_h', T.Properties.VariableNames)
        BCS.SDRidgePeriod = double(T.IQR_RidgePeriod_h);
    else
        BCS.SDRidgePeriod = nan(n, 1);
    end
    BCS.MeanRidgePower_log10 = double(T.MeanRidgePower_log10);
    BCS.SDRidgePower_log10 = double(T.SDRidgePower_log10);
end
