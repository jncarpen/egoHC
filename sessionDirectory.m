function [animal_trackfilelist, animal_fmEvents, animal_behaviourlist, animal_posfilelist, animal_tfilelist, animal_EEGfilelist] = sessionDirectory(recFolderPath)
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

    animal_tfilelist = cell(1,length(dirList));
    animal_EEGfilelist = cell(1,length(dirList));
    animal_posfilelist = cell(1,length(dirList));
    animal_trackfilelist = cell(1,length(dirList));
    animal_fmEvents = cell(1,length(dirList));
    animal_behaviourlist = cell(1, length(dirList));


for session= 1:length(dirList)
    folderPath = char(dirList(session));
    
    % grab files with specific extensions
    fileList_t64 = dir(fullfile(folderPath, '*.t64'));
    tfilelist = cell(1,length(fileList_t64));

    fileList_pos = dir(fullfile(folderPath, '*.pos'));
    fileList_tracker = dir(fullfile(folderPath, '*_tracker.mat'));
    fileList_fmEvents = dir(fullfile(folderPath, 'fmEvents.mat'));
    fileList_behaviour = dir(fullfile(folderPath, 'behaviourTimes.mat'));
    
%     % if no file is found...
%     if isempty(fileList_tracker)
%         fileList_tracker = [];
%     end
%     
%     if isempty(fileList_fmEvents)
%         fileList_fmEvents = [];
%     end
    
    % find a less gross way to do this 
    pathEEG = dir(fullfile(folderPath, '*.eeg'));
    pathEEG2 = dir(fullfile(folderPath, '*.eeg2'));
    pathEEG3 = dir(fullfile(folderPath, '*.eeg3'));
    pathEEG4 = dir(fullfile(folderPath, '*.eeg4'));
    fileList_EEG = [pathEEG, pathEEG2, pathEEG3, pathEEG4];
    EEGfilelist = cell(1,length(fileList_EEG));
        
    for i=1:length(fileList_t64)
        tfilelist{i} = append(fileList_t64(i).folder, '\', fileList_t64(i).name);
    end
    
    for i=1:length(fileList_EEG)
        EEGfilelist{i} = append(fileList_EEG(i).folder, '\', fileList_EEG(i).name); 
    end
    
    for i=1:length(fileList_pos)
        posfilelist{i} = append(fileList_pos(i).folder, '\', fileList_pos(i).name); 
    end
    
     trackerfilelist{1}=[];
    for i=1:length(fileList_tracker)
        if ~isempty(fileList_tracker(i))
            trackerfilelist{i} = append(fileList_tracker(i).folder, '\', fileList_tracker(i).name); 
        else
            trackerfilelist{i} = 'EMPTY';
        end
    end
    
     behaviorfilelist{1}=[];
    for i=1:length(fileList_behaviour)
        if ~isempty(fileList_behaviour(i))
            behaviorfilelist{i} = append(fileList_behaviour(i).folder, '\', fileList_behaviour(i).name); 
        else
            behaviorfilelist{i} = 'OFS';
        end
    end
    
    fmEventsfilelist{1}="OFS";
    for i=1:length(fileList_fmEvents)
        if ~isempty(fileList_fmEvents(i))
            fmEventsfilelist{i} = append(fileList_fmEvents(i).folder, '\', fileList_fmEvents(i).name); 
        else
            fmEventsfilelist{i} = 'OFS';
        end
    end
    
 animal_tfilelist{session} = tfilelist;
 animal_EEGfilelist{session} = EEGfilelist;
 animal_posfilelist{session} = posfilelist;
 animal_trackfilelist{session} = trackerfilelist;
 animal_fmEvents{session} = fmEventsfilelist;
 animal_behaviourlist{session} = behaviorfilelist;
 
 clear fileList_fmEvents fileList_tracker
 
end

end

