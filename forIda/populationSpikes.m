% folder that you want to pull matfiles from
folderPath =  'C:\Users\17145\Documents\github_local\egoHC\dManagement';

% select matfiles to pull out of the folder
matfilelist = dir(fullfile(folderPath, '*.m'));

% make a list of the paths to these matfiles
for i=1:length(matfilelist)
    pathnames{i} = append(matfilelist(i).folder, '\', matfilelist(i).name);
end

spikearray = cell(1,length(pathnames));
for celnum = 1:length(pathnames)
    S = load(pathnames{1,celnum});
    spikearray{1,celnum} = S.cellTS;
    lengthVec(celnum) = length(S.cellTS);
end