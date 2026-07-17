function E = extended_ll_projected_event_definitions(refPhotoperiodH)
%EXTENDED_LL_PROJECTED_EVENT_DEFINITIONS Projected DL/LD/mid anchors under LL.
%
%   E = extended_ll_projected_event_definitions(refPhotoperiodH)
%
%   Default reference photoperiod is L22 (22 h light). Event ZTs:
%     ProjectedDL_LL      ZT 0
%     ProjectedLD_LL      ZT refPhotoperiodH (lights-off of L22 = ZT22)
%     ProjectedMidLight   mid of projected light
%     ProjectedMidDark    mid of projected dark (ZT22-24 under L22)
%
%   Extracted from Kent D v4 LL projected aftereffect analysis.
    refDarkH = 24 - refPhotoperiodH;
    types = ["ProjectedDL_LL"; "ProjectedLD_LL"; "ProjectedMidLight_LL"; "ProjectedMidDark_LL"];
    baseTypes = ["DL"; "LD"; "MidLight"; "MidDark"];
    zts = [0; refPhotoperiodH; refPhotoperiodH/2; mod(refPhotoperiodH + refDarkH/2, 24)];
    E = table(types, baseTypes, zts, 'VariableNames', {'TransitionType','BaseTransitionType','TransitionZT_h'});
end

