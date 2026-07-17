function tf = extended_period_gate_same_key(T, r)
%EXTENDED_PERIOD_GATE_SAME_KEY Match File/SignalID/Condition/Photoperiod/Phase/Band.

    tf = true(height(T), 1);
    tf = tf & strcmpi(string(T.File), string(r.File));
    tf = tf & strcmpi(string(T.SignalID), string(r.SignalID));
    tf = tf & strcmpi(string(T.ConditionParsed), string(r.ConditionParsed));
    tf = tf & abs(T.Photoperiod_h - r.Photoperiod_h) < 1e-9;
    tf = tf & strcmpi(string(T.Phase), string(r.Phase));
    tf = tf & strcmpi(string(T.BandName), string(r.BandName));
end
