function extended_period_gate_delete_if_exists(p)
%EXTENDED_PERIOD_GATE_DELETE_IF_EXISTS Delete file if present.
    if isfile(p), delete(p); end
end
