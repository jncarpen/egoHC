function [distHome] = objDist(pos, events)
%OBJDIST Summary of this function goes here
%   INPUT:
%   pos:        matrix in the form [t x y]
%   events:     struct w/ fields, where N is # of events
%                   (1) events.times (Nx2 double)
%                   (2) events.type (Nx1 cell)
%                   (3) events. locations (Nx1 double)
%                   (4) events.dataPath (string)


    %% compute distance from home well

    % Parse position vector
    t = pos(:,1); % time
    x = pos(:,2); % posx
    y = pos(:,3); % posy

    % find home well location
    homeLocIdx = find(events.type == "DRINKING_HOME", 1, 'first');
    homeLoc = events.locations(homeLocIdx);

    % Convert FM well location index to [x,y] coordinates \\ 
        % Need to find out where to get this information.
    homeXY(1,1:2) = [,];

    % Calculate distance from LED to *home well* for each timestamp.
        % Change out ObjCoord_2 for whatever the well location is
    for cordInd = 1:length(x)
        DisValue = sqrt((homeXY(1,1)- x(cordInd))^2+ (homeXY(1,2)- y(cordInd))^2); % distance formula
        distHome(cordInd) = DisValue;
    end

end