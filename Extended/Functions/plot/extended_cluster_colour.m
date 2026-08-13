function rgb = extended_cluster_colour(pal, varargin)
%EXTENDED_CLUSTER_COLOUR Tol Bright RGB for locked cluster mean traces / violins.
%
%   rgb = extended_cluster_colour(pal, 'BandName', "UR_3_6", 'ClusterRank', 2)
%   rgb = extended_cluster_colour(pal, 'ClusterID', cid, 'ClusterSummary', CS)
%
%   Locked keys (BandName|ClusterRank):
%     UR_1_3|1  — Cyan   (#66CCEE)  C01
%     UR_3_6|1  — Purple (#AA3377)  C01
%     UR_3_6|2  — Red    (#EE6677)  C02
%
%   Sex colours (green/yellow) remain reserved for F|M splits only.
%   Individual mouse traces stay light grey ([0.70–0.75]).

    p = inputParser;
    addParameter(p, 'BandName', "", @(x) true);
    addParameter(p, 'ClusterRank', NaN, @(x) true);
    addParameter(p, 'ClusterID', "", @(x) true);
    addParameter(p, 'ClusterSummary', table(), @(x) istable(x));
    parse(p, varargin{:});

    bandName = string(p.Results.BandName);
    rank = double(p.Results.ClusterRank);
    cid = string(p.Results.ClusterID);
    CS = p.Results.ClusterSummary;

    if strlength(cid) > 0 && ~isempty(CS) && ismember('ClusterID', CS.Properties.VariableNames)
        row = CS(string(CS.ClusterID) == cid, :);
        if ~isempty(row)
            if ismember('BandName', row.Properties.VariableNames)
                bandName = string(row.BandName(1));
            end
            if ismember('ClusterRank', row.Properties.VariableNames)
                rank = double(row.ClusterRank(1));
            end
        end
    end

    key = extended_cluster_key_(bandName, rank);
    if isfield(pal, 'cluster') && isa(pal.cluster, 'containers.Map') && isKey(pal.cluster, key)
        rgb = pal.cluster(key);
        return;
    end

    if isfield(pal, 'band') && isa(pal.band, 'containers.Map')
        bkey = char(string(bandName));
        if isKey(pal.band, bkey)
            rgb = pal.band(bkey);
            return;
        end
    end

    if isfield(pal, 'base') && ~isempty(pal.base)
        rgb = pal.base(7, :);
    else
        rgb = [0.73 0.73 0.73];
    end
end

function key = extended_cluster_key_(bandName, rank)
    bn = char(string(bandName));
    if ~isfinite(rank) || rank <= 0
        key = bn;
    else
        key = sprintf('%s|%d', bn, round(rank));
    end
end
