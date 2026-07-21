function paths = extended_script8_resolve_paths(cohortRoot)
%EXTENDED_SCRIPT8_RESOLVE_PATHS Resolve Script 1-7 inputs from cohort root (e.g. C57_LP).

    cohortRoot = char(string(cohortRoot));
    if ~isfolder(cohortRoot)
        error('extended_script8_resolve_paths:BadRoot', 'Not a folder: %s', cohortRoot);
    end

    [~, cohortTag] = fileparts(cohortRoot);
    paths = struct();
    paths.cohortRoot = cohortRoot;
    paths.cohortTag = cohortTag;

    paths.script3Root = find_first_dir_(cohortRoot, '03_AcrossPhotoperiod_*');
    paths.handoffRoot = fullfile(cohortRoot, 'Handoff');
    paths.extendedHandoff = fullfile(cohortRoot, 'ExtendedHandoff', 'AcrossPhotoperiod_Input');

    paths.scalogramRawF = fullfile(paths.script3Root, 'Figures', 'Stitched_Scalograms', ...
        ['Stitched_Scalogram_Average_F_' cohortTag '.jpeg']);
    paths.scalogramRawM = fullfile(paths.script3Root, 'Figures', 'Stitched_Scalograms', ...
        ['Stitched_Scalogram_Average_M_' cohortTag '.jpeg']);
    paths.scalogramHSubF = fullfile(paths.script3Root, 'Figures', 'Stitched_Scalograms_HSub', ...
        ['Stitched_HSub_Removed_Residual_SEL_P360_F_' cohortTag '.jpeg']);
    paths.scalogramHSubM = fullfile(paths.script3Root, 'Figures', 'Stitched_Scalograms_HSub', ...
        ['Stitched_HSub_Removed_Residual_SEL_P360_M_' cohortTag '.jpeg']);

    paths.gateXlsx = fullfile(paths.extendedHandoff, 'RawVsSelectiveHSub_PeriodValidation', ...
        'HSubSupported_PeriodMap.xlsx');
    paths.gateRetentionFig = fullfile(paths.extendedHandoff, ...
        'RawVsSelectiveHSub_PeriodValidation', 'Figures', 'Retention_ByBand.png');

    paths.resyncXlsx = fullfile(paths.extendedHandoff, 'Ultradian_RidgePhase_Resync', ...
        'Ultradian_RidgePhase_Resync_Output.xlsx');
    paths.resyncMat = fullfile(paths.extendedHandoff, 'Ultradian_RidgePhase_Resync', ...
        'Ultradian_RidgePhase_Resync_Output.mat');
    paths.resyncGradientXlsx = fullfile(paths.extendedHandoff, 'Ultradian_RidgePhase_Resync', ...
        'TransitionEffect_vs_Photoperiod', 'TransitionEffect_vs_Photoperiod.xlsx');

    paths.lmeDescriptiveXlsx = fullfile(paths.extendedHandoff, 'AcrossPhotoperiod_LME', 'Tables', ...
        'AcrossPhotoperiod_Outputs.xlsx');
    paths.lmeInferenceXlsx = fullfile(paths.extendedHandoff, 'AcrossPhotoperiod_LME', 'Tables', ...
        'AcrossPhotoperiod_LME_Outputs.xlsx');

    paths.profilesRoot = pick_profiles_root_(cohortRoot);
    paths.profilesXlsx = fullfile(paths.profilesRoot, 'Tables', ...
        'SelectedValidatedUR_PublicationProfiles_Output.xlsx');
    if ~isfile(paths.profilesXlsx)
        alt = fullfile(paths.profilesRoot, 'Tables', ...
            'SelectedValidatedUR_DevelopmentProfiles_Output.xlsx');
        if isfile(alt)
            paths.profilesXlsx = alt;
        end
    end

    paths.missing = strings(0, 1);
    req = {
        'script3Root', paths.script3Root;
        'scalogramRawF', paths.scalogramRawF;
        'scalogramRawM', paths.scalogramRawM;
        'scalogramHSubF', paths.scalogramHSubF;
        'scalogramHSubM', paths.scalogramHSubM;
        'gateXlsx', paths.gateXlsx;
        'resyncXlsx', paths.resyncXlsx;
        'lmeDescriptiveXlsx', paths.lmeDescriptiveXlsx;
        'lmeInferenceXlsx', paths.lmeInferenceXlsx;
        'profilesXlsx', paths.profilesXlsx;
        };
    for i = 1:size(req, 1)
        p = req{i, 2};
        if isempty(p) || (~isfolder(p) && ~isfile(p))
            paths.missing(end + 1, 1) = string(req{i, 1}); %#ok<AGROW>
        end
    end
end

function d = find_first_dir_(root, pattern)
    d = '';
    hits = dir(fullfile(root, pattern));
    hits = hits([hits.isdir]);
    if ~isempty(hits)
        d = fullfile(hits(1).folder, hits(1).name);
    end
end

function root = pick_profiles_root_(cohortRoot)
    parent = fullfile(cohortRoot, 'ExtendedHandoff');
    pub = fullfile(parent, 'SelectedValidatedUR_PublicationProfiles');
    dev = fullfile(parent, 'SelectedValidatedUR_DevelopmentProfiles');
    if isfolder(pub)
        root = pub;
    elseif isfolder(dev)
        root = dev;
    else
        root = pub;
    end
end
