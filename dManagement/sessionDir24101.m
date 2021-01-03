function [outputArg1,outputArg2] = sessionDir24101(inputArg1,inputArg2)
%SESSIONDIR24101 * Session directory for Rat #24101
%   This should be a temporary script, but accounts for the organizational
%   structure of recordings in animal 24101.
%   
%   INPUTS
%   OUTPUTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recDir = dir(recFolderPath);

count = 1;
for dirIndex = 1:length(recDir)
    L = strlength(recDir(dirIndex).name);
    % check that folder is named as a session
    
    if L == 10 
        % squash strings together to get correct filepath
        dirList{count} = append(recDir(dirIndex).folder, '\', recDir(dirIndex).name);
        count = count + 1;
    end
end

dirList = dirList(~cellfun(@isempty, dirList)); % clear empty cells

%% get correct order for dirList (sorted by date)

MDYS = cell(length(dirList), 1);

for fileNum = 1:length(dirList)
    splitName = split(dirList{1,fileNum},'\');
    serial = splitName{end,1}; % grab serial #
    MDYS{fileNum, 1} = strcat(serial(3:4), serial(1:2), serial(5:8), serial(9:10));
end

% sort by mo-day-yr-session
[~, I] = sort(MDYS); % grab order indices
dirList = dirList(I); 

%% 







end

