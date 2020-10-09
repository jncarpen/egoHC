% Fix sessType (as generated in the generateDataset.m script)
% This can be deleted once I fix the other script
% July 23, 2020. Jordan Carpenter 

trialType = cell(1, length(sessType_A)); % initalize new cell array

for sess = 1:length(sessType_A)
    if ~isempty(sessType_A{1,sess})
        if string(sessType_A{1,sess}) == "Fm trial" || string(sessType_A{1,sess}) == "Fm" || string(sessType_A{1,sess}) == "fmc" 
            trialType{1,sess} = 'FM';
        elseif string(sessType_A{1,sess}) == "Of"
            trialType{1,sess} = 'OF';
        end
    end
end
