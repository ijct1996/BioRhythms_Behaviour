function extended_period_gate_write(outXLSX, outMAT, validationMap)
%EXTENDED_PERIOD_GATE_WRITE Write HSubSupported_PeriodMap.xlsx and .mat.

    extended_period_gate_delete_if_exists(outXLSX);

    Settings = validationMap.Settings;
    LoadSummary = validationMap.LoadSummary;
    allCandidates = validationMap.AllPeriodCandidates;
    Matched = validationMap.Matched_Periods_All;
    CarryForward = validationMap.CarryForward_Periods;
    RawOnly = validationMap.RawOnly_NotCarriedForward;
    HSubOnly = validationMap.HSubOnly_ResidualFeatures;
    FLSens = validationMap.FullLadder_Sensitivity;
    RetBand = validationMap.Retention_ByBand;
    RetPhoto = validationMap.Retention_ByPhotoperiod;
    RetPhotoBand = validationMap.Retention_ByPhotoperiodBand;
    QCFlags = validationMap.QC_Flags;

    readme = {
        'Raw vs Selective-HSub ultradian period validation.'
        'Primary rule: Raw ultradian candidates are carried forward only when matched to a SEL_P360 HSub candidate within the configured period tolerance.'
        'Raw remains the biological signal. HSub is used only as validation.'
        'Full Ladder outputs are retained as sensitivity/stress-test flags, not as the primary validation rule.'
        'UR_12_18 is flagged as harmonic-sensitive because it is close to the first harmonic of a 24 h rhythm.'
        'CarryForward_Periods is the table to use in downstream validated-Raw analyses.'
        };
    writecell(readme(:), outXLSX, 'Sheet', 'README');
    writetable(Settings, outXLSX, 'Sheet', 'Settings');
    writetable(LoadSummary, outXLSX, 'Sheet', 'LoadSummary');
    writetable(allCandidates, outXLSX, 'Sheet', 'AllPeriodCandidates');
    writetable(Matched, outXLSX, 'Sheet', 'Matched_Periods_All');
    writetable(CarryForward, outXLSX, 'Sheet', 'CarryForward_Periods');
    writetable(RawOnly, outXLSX, 'Sheet', 'RawOnly_NotCarriedForward');
    writetable(HSubOnly, outXLSX, 'Sheet', 'HSubOnly_ResidualFeatures');
    writetable(FLSens, outXLSX, 'Sheet', 'FullLadder_Sensitivity');
    writetable(RetBand, outXLSX, 'Sheet', 'Retention_ByBand');
    writetable(RetPhoto, outXLSX, 'Sheet', 'Retention_ByPhotoperiod');
    writetable(RetPhotoBand, outXLSX, 'Sheet', 'Retention_ByPhotoBand');
    writetable(QCFlags, outXLSX, 'Sheet', 'QC_Flags');

    save(outMAT, 'validationMap', '-v7.3');
end
