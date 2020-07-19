function [unitInfo] = pullLabNotes(LabNotes_flnm, totalTetrodes)
%PULLLABNOTES 
%   INPUT:
%   LabNotes_flnm:  .xlsx sheet filepath for a *SINGLE* animal
%   totalTetrodes:  Total # of tetrodes that animal had implanted.

%   OUTPUT:
%   unitInfo:       table w/ 6 variables [CutFile, TypeOfSession, Tetrode, Cluster, UniqueID, ID]

%   jCarpenter

tetNum = ["TT01", "TT02", "TT03", "TT04", "TT05", "TT06", "TT07", "TT08"]; % add as many as needed...
% Tet = cell(1,totalTetrodes);

if totalTetrodes > length(tetNum)
    disp("Error: Add options for additional tetrodes in getUniqueID.m")
end

unitInfo = []; % initialize empty 'table'

for tetSheet = 1:totalTetrodes
tetInfo = readtable(LabNotes_flnm, 'Sheet', tetNum(tetSheet));

    cRows = [];
    for row = 1:height(tetInfo)
        if isnan(tetInfo.Cluster(row))
            cRows(row) = row;
        end
    end
    cRows = nonzeros(cRows); 
    tetInfo(cRows,:) = []; % delete 'c' rows
    
    % Columns to keep 
    tetInfoTite = table(tetInfo.CutFile, tetInfo.TypeOfSession, tetInfo.Tetrode, tetInfo.Cluster, tetInfo.UniqueID, tetInfo.ID);
    tetInfoTite.Properties.VariableNames = {'CutFile', 'TypeOfSession', 'Tetrode', 'Cluster', 'UniqueID', 'ID'};

    for row = 1:height(tetInfoTite)
        cutFl = tetInfoTite.CutFile{row,1};
        cutFl = extractAfter(cutFl,"recordings\"); % I want to remove the .clusters at the end but not every string has it
        if strfind(cutFl,".clusters") == 13
            cutFl = extractBefore(cutFl,".clusters");
        end
        
        if strfind(cutFl, "_") == 11 % in position 11?
            cutFl = extractBefore(cutFl, "_");
        end
        tetInfoTite.CutFile{row,1} = cutFl;
    end
    
    unitInfo = [unitInfo; tetInfoTite]; % concat tables from multiple tetrodes together
end

end

