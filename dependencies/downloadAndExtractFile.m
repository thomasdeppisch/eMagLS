function downloadAndExtractFile(file, url)
    % Download a file from specified URL to the designated target location.
    % 
    % If the provided file in the URL as a ZIP archive, then a file with a
    % name identical to the target will be extracted.
    %
    % Hannes Helmholz, 2023

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
