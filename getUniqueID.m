function [uniqueID_A, neuronID_A, sessType_A] = getUniqueID(unitInfo, animal_tfilelist)
%GETUNIQUEID Find uniqueID for each cell for an animal.
%   INPUT:
%   unitInfo:               can be generated with function pullLabNotes.m *
%   animal_tfilelist:       filelist of .t64 file locations

%   OUTPUT:
%   uniqueID:       unique cell ID (should be in the same format as S
%                   [cellspikes])
%   sessType:       session type of each cell

% jCarpenter

%% FUNCTION

uniqueID_A = cell(1,length(animal_tfilelist));
neuronID_A = cell(1,length(animal_tfilelist));
sessType_A = cell(1,length(animal_tfilelist));

for sessNum = 1:length(animal_tfilelist)

t64List = animal_tfilelist{1,sessNum};
numCells = length(t64List); % list of t64 file locations for a single session **
t64info = cell(1, numCells);
uniqueID = cell(1,numCells);
neuronID = cell(1,numCells);
sessType = cell(1,numCells);

for neuron = 1:numCells
    t64str = t64List{1,neuron};
    splitPath = split(t64str,'\');
    t64name = splitPath{end,1};
    t64name = extractBefore(t64name,".t64");
    t64name = split(t64name,'_'); % [day-mo-yr-trial# (serial), tet#, clust#]
    t64name{3,1} = regexprep(t64name{3,1},'^0*',''); % pull out leading zeros 
    t64info{1,neuron} = t64name;
end

for neuron = 1:length(t64info)
    serial_t64 = convertCharsToStrings(t64info{1,neuron}{1,1});
    tetrode_t64 = str2double(convertCharsToStrings(t64info{1,neuron}{2,1}));
    cluster_t64 = str2double(convertCharsToStrings(t64info{1,neuron}{3,1}));
    
    for row = 1:height(unitInfo)
        serial_unit = convertCharsToStrings(unitInfo.CutFile{row,1});
        tetrode_unit = convertCharsToStrings(unitInfo.Tetrode(row,1));
        cluster_unit = convertCharsToStrings(unitInfo.Cluster(row,1));
        unique = unitInfo.UniqueID(row,1);
        ID = unitInfo.ID(row,1);
        type = unitInfo.TypeOfSession{row,1};
        
        if serial_t64 == serial_unit && tetrode_t64 == tetrode_unit && cluster_t64 == cluster_unit
            uniqueID{1,neuron} = unique;
            neuronID{1,neuron} = ID;            
        end
        
        if serial_t64 == serial_unit
            sessType = type;
        end
    end
end


uniqueID_A{1,sessNum} = uniqueID;
neuronID_A{1,sessNum} = neuronID;
sessType_A{1,sessNum} = sessType;

end
end % end function
