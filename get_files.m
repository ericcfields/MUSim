%Get files in a given directory with a given extension
%
%INPUTS
% directory  - directory containing files
% extension  - file extentsion (e.g., '.txt') 
%              {default: return all files (but not directories)}
%
%OUTPUT
% files      - a cell array of files
%
%Author: Eric Fields
%Version Date: 12 March 2019

function files = get_files(directory, extension)
    
    %Get cell array of files
    files = dir(directory);
    files = files(~[files.isdir]); %remove directories
    files = {files.name}; %get cell array of file names
    
    %Find subset with given extension if requested
    if nargin > 1 && ~isempty(extension)
        idx = cellfun(@(x) char_endswith(x, extension, false), files);
        files = files(idx);
    end
    
end

function TF = char_endswith(text, pattern, match_case)
%Return boolean indicating whether a string (text) ends with another 
%string (pattern)
    
    %Get ending of appropriate length
    if length(text) >= length(pattern)
        text_end = text((end-length(pattern)+1):end);
    else
        TF = false;
        return;
    end
    
    %Default to exact match
    if nargin < 3
        match_case = true;
    end
    
    %Compare
    if match_case
        TF = strcmp(text_end, pattern);
    else
        TF = strcmpi(text_end, pattern);
    end

end
