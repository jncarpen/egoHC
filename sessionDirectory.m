function animal_tfilelist = sessionDirectory(recFolderPath)
%SESSIONDIRECTORY:

    % INPUT:
    % recFolderPath: master folder for a single animal's recording sessions should be read in as
    % (i.e. 'D:\Data\24116\Axona recordings'). 
    
    % OUTPUT:
    % animal_tfilelist: 1 x S cell array, where S is the number of sessions
    % recorded for the animal. Each cell contains a 1 x T cell array where
    % T is the number of .t64 files found in each session folder. Each cell
    % contains information about the file path of each .t64 file to be read
    % into AD Redish's 'LoadSpikes.m' script to pull out spike trains.

    % Jo Carpenter, June 23, 2020

recDir = dir(recFolderPath);
dirList = cell(1, length(recDir));
animal_tfilelist = cell(1,length(recDir));

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

for session= 1:length(dirList)
    folderPath = char(dirList(session));
    % grab files with .t64 extension
    fileList = dir(fullfile(folderPath, '*.t64'));
    tfilelist = cell(1,length(fileList));
    
    for i=1:length(fileList)
        tfilelist{i} = append(fileList(i).folder, '\', fileList(i).name);
    end
    
 animal_tfilelist{session} = tfilelist;
 
end
end

