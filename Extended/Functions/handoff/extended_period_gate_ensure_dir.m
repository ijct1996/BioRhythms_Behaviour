function extended_period_gate_ensure_dir(p)
%EXTENDED_PERIOD_GATE_ENSURE_DIR Create folder if missing.
    if ~isfolder(p), mkdir(p); end
end
