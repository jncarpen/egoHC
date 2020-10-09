function defineWells(trialType, pos)
%DEFINEWELLS Summary of this function goes here
%   INPUTS
%   trialType:      1 x S cell array, where S is the number of sessions.
%   pos:            [t x1 y1 x2 y2]
% 
%   J. Carpenter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for sess = 1:10%length(trialType)
    if trialType{1,sess} == 'FM'
        plot(pos{1,sess}(:,2),pos{1,sess}(:,3), 'k')
        hold on
        % pause
    end
end

