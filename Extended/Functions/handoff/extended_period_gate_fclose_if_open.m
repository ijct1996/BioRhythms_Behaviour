function extended_period_gate_fclose_if_open(fid)
%EXTENDED_PERIOD_GATE_FCLOSE_IF_OPEN Close file id if still open.
    if ~isempty(fid) && fid > 0
        fclose(fid);
    end
end
