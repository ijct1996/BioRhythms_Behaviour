function id = extended_ridge_make_candidate_id(fileStem, signalID, sourceTag, hsubModeTag, phaseTag, bandName, rank)
    raw = sprintf('%s__%s__%s__%s__%s__%s__R%d', ...
        char(string(fileStem)), char(string(signalID)), char(string(sourceTag)), ...
        char(string(hsubModeTag)), char(string(phaseTag)), char(string(bandName)), rank);
    id = extended_ridge_sanitise_filename(raw);
end

function s = extended_ridge_sanitise_filename(strIn)
    s = regexprep(string(strIn), '[^\w\-]', '_');
    s = char(s);
    if isempty(s), s = 'unnamed'; end
end
