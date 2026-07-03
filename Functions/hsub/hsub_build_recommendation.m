function rec = hsub_build_recommendation(mouseOut, opts)
%HSUB_BUILD_RECOMMENDATION Per-mouse residual recommendation (default SEL_P360).

    cfg = core_defaults();
    if nargin < 2 || isempty(opts)
        opts = hsub_get_opts();
    end
    opts.defaultResidual = cfg.hsub.defaultResidual;

    rec = struct('code', 'PASS', 'mode', 'Selective', 'minKey', 'Min360', ...
        'anchorOK', false, 'reason', 'Anchor rejected — pass-through raw.');

    if ~mouseOut.anchorOK
        return;
    end

    rec.anchorOK = true;
    chosenMinKey = 'Min360';
    if isfield(mouseOut, 'crossQC') && isfield(mouseOut.crossQC, 'Min60')
        nRMS = mouseOut.crossQC.Min60.nRMS;
        if isfinite(nRMS) && nRMS >= opts.nRMS_high
            chosenMinKey = 'Min60';
        end
    end

    mode = 'Selective';
    if isfield(mouseOut.crossQC, chosenMinKey)
        qc = mouseOut.crossQC.(chosenMinKey);
        if qc.nRMS >= opts.nRMS_high || qc.dVar >= opts.dVar_high
            mode = 'FullLadder';
        end
    end

    rec.minKey = chosenMinKey;
    rec.mode = mode;

    switch mode
        case 'FullLadder'
            rec.code = sprintf('FULL_%s', strrep(chosenMinKey, 'Min', 'P'));
        otherwise
            rec.code = sprintf('SEL_%s', strrep(chosenMinKey, 'Min', 'P'));
    end

    rec.reason = sprintf('%s: mode=%s, %s.', rec.code, mode, chosenMinKey);
end
