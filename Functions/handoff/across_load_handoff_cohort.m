function out = across_load_handoff_cohort(handoffDir)
%ACROSS_LOAD_HANDOFF_COHORT Load all CoreSummary MAT files in Handoff/.

    mats = dir(fullfile(handoffDir, 'CoreSummary__*.mat'));
    if isempty(mats)
        error('across_load_handoff_cohort:NoFiles', 'No CoreSummary__*.mat in %s', handoffDir);
    end

    out = struct('fileStem', {}, 'lightHours', {}, 'periodTable', {}, 'bandPower', {}, ...
        'meta', {}, 'groups', {}, 'paths', {}, 'summaryPath', {});
    for i = 1:numel(mats)
        matPath = fullfile(mats(i).folder, mats(i).name);
        S = load(matPath);
        out(i).fileStem = S.fileStem;
        out(i).lightHours = S.meta.lightDurationHours;
        out(i).periodTable = S.periodTable;
        out(i).bandPower = S.bandPower;
        out(i).meta = S.meta;
        out(i).hsub = S.hsub;
        out(i).groups = S.groups;
        out(i).paths = S.paths;
        out(i).summaryPath = matPath;
    end

    [~, ord] = sort([out.lightHours]);
    out = out(ord);
end
