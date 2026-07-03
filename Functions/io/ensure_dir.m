function ensure_dir(dirPath)
%ENSURE_DIR Create directory if missing.
    if ~isfolder(dirPath)
        mkdir(dirPath);
    end
end
