function [recDir, animal_tfilelist, animal_EEGfilelist] = sessionDirectory(recFolderPath)
%SESSIONDIRECTORY:

    % INPUT:
    % recFolderPath: master folder for a single animal's recording sessions should be read in as
    % (i.e. 'D:\Data\24116\Axona recordings'). 
    
    % OUTPUT:
    % (1) animal_tfilelist: 1 x S cell array, where S is the number of sessions
    % recorded for the animal. Each cell contains a 1 x T cell array where
    % T is the number of .t64 files found in each session folder. Each cell
    % contains information about the file path of each .t64 file to be read
    % into AD Redish's 'LoadSpikes.m' script to pull out spike trains.
    % Sessions without .t64 files will maintain their position in
    % animal_tfilelist but will be empty.
    % (2) animal_EEGfilelist: paths to all four EEG file locations for each sesssion
    % (Axona only)- will be adapted for Neuralynx.
    % (3) recDir: Directory of recording folder paths.

    % Jordan Carpenter, June 23, 2020

    recDir = dir(recFolderPath);
    dirList = cell(1, length(recDir));
    animal_tfilelist = cell(1,length(recDir));
    animal_EEGfilelist = cell(1, length(recDir));

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

% % remove junk rows from recDir struct
% dirName = [recDir.name]; % Get all the names
% g = length(dirName) == 10; % find names with correct format
% recDir = recDir(g); % Select only correct names, delete rest.

dirList = dirList(~cellfun(@isempty, dirList)); % clear empty cells

for session= 1:length(dirList)
    folderPath = char(dirList(session));
    
    % grab files with .t64 or .eeg extension
    fileList_t64 = dir(fullfile(folderPath, '*.t64'));
    
    % find a less gross way to do this
    pathEEG = dir(fullfile(folderPath, '*.eeg'));
    pathEEG2 = dir(fullfile(folderPath, '*.eeg2'));
    pathEEG3 = dir(fullfile(folderPath, '*.eeg3'));
    pathEEG4 = dir(fullfile(folderPath, '*.eeg4'));
    fileList_EEG = [pathEEG, pathEEG2, pathEEG3, pathEEG4];
    
    tfilelist = cell(1,length(fileList_t64));
    EEGfilelist = cell(1,length(fileList_EEG));
    
    for i=1:length(fileList_t64)
        tfilelist{i} = append(fileList_t64(i).folder, '\', fileList_t64(i).name);
    end
    
    for i=1:length(fileList_EEG)
        EEGfilelist{i} = append(fileList_EEG(i).folder, '\', fileList_EEG(i).name); 
    end
    
 animal_tfilelist{session} = tfilelist;
 animal_EEGfilelist{session} = EEGfilelist;
 
end
end

