function [idLookup] = cell_index(UniqueID, unitInfo)
%CELL_INDEX Make a list of indices for each unique unit.
%   
%   USAGE               [idLookup] = cell_index(UniqueID, unitInfo);
%
%   INPUTS
%   UniqueID:           1xS cell array, where S is the number of sessions for a
%                       given animal. Each of the S cell arrays contains
%                       another cell array of length 1xU, where U is the
%                       number of units recorded during each session. This
%                       input can be generated with the getUniqueID.m
%                       function.
%
%   unitInfo:           Jx6 table, where J is the number of non-unique unit
%                       identifiers for a given animal. The important
%                       information for this function is the
%                       unitInfo.UniqueID column, which will be used to
%                       compute the highest-valued unique identifier.
%   OUTPUTS
%   idLookup:           Ux2 cell array, where each row corresponds to a
%                       unique cell ID. The second column is a list of Sess/Unit
%                       pairs: [Unique ID, [Sess, Unit]].
%   Jordan Carpenter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalSessions = length(UniqueID);
maxID = nanmax(unitInfo.UniqueID);
idLookup = cell(maxID, 2);

for ID = 1:maxID
    actvCells = [];
    for sess = 1:totalSessions
        for unit = 1:length(UniqueID{1,sess})
            if UniqueID{1,sess}{1,unit} == ID
                actvCells = [actvCells; [sess, unit]];
            end
        end
    end
    idLookup{ID,1} = ID;
    idLookup{ID,2} = actvCells;   
end
end

