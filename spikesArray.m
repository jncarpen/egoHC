function [cellSpikes, unitID] = spikesArray(animal_tfilelist)
%SPIKESARRAY Summary of this function goes here
%   Detailed explanation goes here

% get length of sessions
totalSessions = length(animal_tfilelist);

% create cell array, where each cell is a session.
cellSpikes = cell(1,totalSessions);
unitID = cell(1,totalSessions);

% loop through the filelists for each session
for sessionNum = 1:totalSessions
    [S, ID] = LoadSpikes(animal_tfilelist{1,sessionNum});
    cellSpikes{1, sessionNum} = S;
    unitID{1,sessionNum} = ID;
end
end

