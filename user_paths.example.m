function up = user_paths()
%USER_PATHS Machine-local output roots (copy to user_paths.m and edit).
%
%   up = user_paths()
%
%   This file is gitignored. Copy from user_paths.example.m:
%       copyfile('user_paths.example.m', 'user_paths.m')
%
%   Fields:
%       outputRoot  — default parent folder for Core v1 results
%       lastCohort  — last cohort tag used (e.g. 'C57_LP')

    up.outputRoot = 'C:\ResearchData\Results\Core';
    up.lastCohort = '';

    % Optional: remember last input folder per script
    up.lastHSubInput  = '';
    up.lastWaveInput  = '';
    up.lastAcrossInput = '';
end
