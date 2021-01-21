% grab # of sessions
totalSessions = length(animal_tfilelist);
stVec = cell(1,totalSessions); 

% grab all the stuff
for sessionNum = 1:totalSessions
    stVec{1, sessionNum} = LoadSpikes(animal_tfilelist{1,sessionNum});
end

% put spikes in array for each session
STV = cell(1,totalSessions); 
for i = 1:length(stVec)
    spkNow = cell(1,length(stVec{1,i}));
    if ~isempty(stVec{1,i})
        for j = 1:length(stVec{1,i})
            spkNow{1,j} = stVec{1,i}{j}.T;
        end
    else
        spkNow{1,j} = [];
    end
    STV{1,i} = spkNow;
end



