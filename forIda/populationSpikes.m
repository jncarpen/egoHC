% folder that you want to pull matfiles from
folderPath =  'C:\Users\17145\Documents\github_local\egoHC\dManagement';

% select matfiles to pull out of the folder
matfilelist = dir(fullfile(folderPath, '*.mat'));

% make a list of the paths to these matfiles
for i=1:length(matfilelist)
    pathnames{i} = append(matfilelist(i).folder, '\', matfilelist(i).name);
end

% load everything in 
spikearray = cell(1,length(pathnames));
for celnum = 1:length(pathnames)
    S = load(pathnames{1,celnum});
    spikearray{1,celnum} = S.cellTS;
    lengthVec(celnum) = length(S.cellTS);
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
populationSpikes = spikesMat';