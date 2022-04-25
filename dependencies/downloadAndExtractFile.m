function downloadAndExtractFile(file, url)
    [filePath, fileName, fileExt] = fileparts(file);
    fileName = [fileName, fileExt];
    [~, ~] = mkdir(filePath); % ignore warning if directory already exists

    % download file
    [~, urlName, urlExt] = fileparts(url);
    urlName = fullfile(filePath, urlName);
    urlFile = [urlName, urlExt];
    if ~isfile(urlFile)
        fprintf('from %s ... ', url);
        websave(urlFile, url);
    end

    if ~isfile(file) 
        % extract ZIP
        if endsWith(urlFile, 'zip')
            fprintf('unpacking archive "%s" ... ', urlFile);
            zipFiles = unzip(urlFile, urlName);
    
            fprintf('extracting "%s" ... ', fileName);
            zipIndex = find(endsWith(zipFiles, fileName));
            if isempty(zipIndex) || length(zipIndex) > 1
                error('Unknown archive content "%s".', fileName);
            end
            movefile(zipFiles{zipIndex}, filePath);
        else
            error('Unknown file format "%s".', urlExt);
        end
    end

    fprintf('done.\n');
end
