%% for ida 

% folder that you want to pull matfiles from
folderPath =  'C:\Users\17145\Documents\github_local\egoHC\dManagement';

% select matfiles to pull out of the folder
matfilelist = dir(fullfile(folderPath, '*.mat'));

% make a list of the paths to these matfiles
for i=1:length(matfilelist)
    pathnames{i} = append(matfilelist(i).folder, '\', matfilelist(i).name);
end

% load everything in and query length of spiketimes
thresh_ceil = 30000;
count = 1; % start a count
for celnum = 1:length(pathnames)
    S = load(pathnames{1,celnum});
    lengthUnitNow = length(S.cellTS);
    if  lengthUnitNow < thresh_ceil
        lengthVec(count) = length(S.cellTS);
        spikearray{1,count} = S.cellTS;
        count = count + 1; % only count units that pass the threshhold
    end
end

% max number of spikes a cell fires during the session
maxspikes = nanmax(lengthVec);

% make a matrix of the correct size
spikesMat = zeros(length(spikearray), maxspikes);

% put spikes into the matrix
for unit = 1:length(spikearray)
    lenNow = lengthVec(unit);
    spikesMat(unit, 1:lenNow) = spikearray{1,unit};
end

% transpose the spikes matrix
popSpikes = spikesMat';

% save the matrix
savepath = 'ida\something';
savename = 'popSpikes_session1.mat';
savestring = strcat(savepath, '\', savename);
save(savestring, popSpikes, '-v7.3')

