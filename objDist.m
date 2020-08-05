function [distHome, spkDist] = objDist(pos, hwLoc, SpikeTimes)
%OBJDIST Summary of this function goes here
%   INPUTS
%   pos:            matrix in the form [t x y]
%   events:         struct w/ fields, where N is # of events
%                       (1) events.times (Nx2 double)
%                       (2) events.type (Nx1 cell)
%                       (3) events. locations (Nx1 double)
%                       (4) events.dataPath (string)
%   hwLoc:          1xS cell array with reward well locations (in FM
%                   coordinates).
%   OUTPUTS
%   distHome:       1xT vector, where T is the number of timestamps in the
%                   pos vector.
%   spkDist:        Distance from home well each time the cell spikes (to
%                   be used in a tuning curve).
%
%   Jordan Carpenter, 2020.
%
% NOTE: Convert distance output to cm **
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% compute distance from home well

    % Parse position vector
    t = pos(:,1); % time
    x = pos(:,2); % posx
    y = pos(:,3); % posy

    % Convert FM well location index to [x,y] coordinates
    switch hwLoc
        case 36
            homeXY(1,1:2) = [372.31, 269.72];
        case 37
            homeXY(1,1:2) = [372.58, 305.75];
    end

    % Calculate distance from LED to *home well* for each timestamp.
        % Change out ObjCoord_2 for whatever the well location is
    for cordInd = 1:length(x)
        DisValue = sqrt((homeXY(1,1)- x(cordInd))^2+ (homeXY(1,2)- y(cordInd))^2); % distance formula
        distHome(cordInd) = DisValue;
    end
    
    %% Get spike distance
    % Compute distance from home well each time the cell spikes.
    
    % find closest timestamp & store values
    idx = knnsearch(distHome', SpikeTimes);
    spkDist = distHome(idx);
    
    if length(SpikeTimes) ~= length(spkDist)
        display("Error: Number of elements in spkDist and SpikeTimes are not the same.")
    end


end