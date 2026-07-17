function T = extended_bh_fdr(T, pVar, alpha)
%EXTENDED_BH_FDR Benjamini-Hochberg FDR within each FDR_Family column.
%
%   T = extended_bh_fdr(T, pVar, alpha)
%
%   Adds Q_BH, BH_Rank, BH_CriticalValue, Significant_BH, Alpha_FDR.
%   Correction is applied separately within each unique FDR_Family.
%   Extracted from Kent D (Ultradian_RidgePhase_Resync_v4_FDR_LLProjected).
    % Benjamini-Hochberg FDR within each FDR_Family.
    if isempty(T) || height(T)==0
        T.Q_BH = nan(height(T),1);
        T.BH_Rank = nan(height(T),1);
        T.BH_CriticalValue = nan(height(T),1);
        T.Significant_BH = false(height(T),1);
        T.Alpha_FDR = repmat(alpha, height(T), 1);
        return;
    end
    if ~ismember(pVar, T.Properties.VariableNames)
        error('P-value column not found for BH/FDR correction: %s', pVar);
    end
    if ~ismember('FDR_Family', T.Properties.VariableNames)
        T.FDR_Family = repmat("Unspecified", height(T), 1);
    else
        T.FDR_Family = string(T.FDR_Family);
    end
    T.Q_BH = nan(height(T),1);
    T.BH_Rank = nan(height(T),1);
    T.BH_CriticalValue = nan(height(T),1);
    T.Significant_BH = false(height(T),1);
    T.Alpha_FDR = repmat(alpha, height(T), 1);

    fams = unique(string(T.FDR_Family), 'stable');
    for f = 1:numel(fams)
        idxFam = find(string(T.FDR_Family)==fams(f));
        p = double(T.(pVar)(idxFam));
        valid = isfinite(p) & p >= 0 & p <= 1;
        if ~any(valid), continue; end
        idxValid = idxFam(valid);
        pValid = p(valid);
        [pSort, ord] = sort(pValid(:), 'ascend');
        m = numel(pSort);
        ranks = (1:m)';
        qSort = pSort .* m ./ ranks;
        qSort = flipud(cummin(flipud(qSort)));
        qSort(qSort > 1) = 1;
        crit = ranks ./ m .* alpha;

        q = nan(m,1); q(ord) = qSort;
        rankOut = nan(m,1); rankOut(ord) = ranks;
        critOut = nan(m,1); critOut(ord) = crit;

        T.Q_BH(idxValid) = q;
        T.BH_Rank(idxValid) = rankOut;
        T.BH_CriticalValue(idxValid) = critOut;
        T.Significant_BH(idxValid) = q <= alpha;
    end
end

