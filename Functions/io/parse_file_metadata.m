function meta = parse_file_metadata(xlsxPath)
%PARSE_FILE_METADATA Parse cohort / photoperiod from Pipeline_Input filename.
%
%   Convention: L{photoperiod}_{genotype}[_LP|_LD|_DD]_02JUL26.xlsx

    [~, stem, ~] = fileparts(xlsxPath);
    meta.filePath = xlsxPath;
    meta.fileStem = stem;
    meta.fileName = [stem '.xlsx'];

    tokens = regexp(stem, '^L(\d+)_(.+)$', 'tokens', 'once');
    if isempty(tokens)
        meta.lightHours = NaN;
        meta.genotypeTag = '';
    else
        meta.lightHours = str2double(tokens{1});
        meta.genotypeTag = tokens{2};
    end

    if contains(stem, '_LP_') || endsWith(stem, '_LP')
        meta.cohort = 'NR2B_LP';
        meta.genotype = 'NR2B';
        meta.protocol = 'LP';
    elseif contains(stem, '_LD_') || endsWith(stem, '_LD')
        meta.cohort = 'NR2B_LD_DD';
        meta.genotype = 'NR2B';
        meta.protocol = 'LD';
    elseif contains(stem, '_DD_') || endsWith(stem, '_DD')
        meta.cohort = 'NR2B_LD_DD';
        meta.genotype = 'NR2B';
        meta.protocol = 'DD';
    elseif contains(stem, '_C57')
        meta.cohort = 'C57_LP';
        meta.genotype = 'C57';
        meta.protocol = 'LP';
    else
        meta.cohort = '';
        meta.genotype = '';
        meta.protocol = '';
    end

    if startsWith(stem, 'L0_')
        meta.lightHours = 0;
    end
end
