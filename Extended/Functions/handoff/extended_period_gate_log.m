function extended_period_gate_log(fid, fmt, varargin)
%EXTENDED_PERIOD_GATE_LOG Print and optionally append a timestamped log line.
    msg = sprintf(fmt, varargin{:});
    fprintf('%s\n', msg);
    if ~isempty(fid) && fid > 0
        fprintf(fid, '[%s] %s\n', datestr(now, 31), msg);
    end
end
