function [list] = AngTrain(STrn, angles)
%ANGTIMES: for making a histogram

list = [];
for numSpk = 1:max(STrn)
    % find indices where this number of spikes occured
    thisManySpikes = find(STrn == numSpk);

    % find values of angular variable at time of spikes
    ang_list = angles(thisManySpikes);

    % replicate this angle (to account for bursting)
    ang_rep = repelem(ang_list, numSpk);
    
    % add to a list
    list = [list; ang_rep'];
end
end

