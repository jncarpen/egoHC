function [homeEvnts, randEvnts] = trigTimes(events)
%PERIEVENT Summary of this function goes here
%   INPUTS
%   events:         struct from fmEvents.mat
%   OUTPUTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = count(events.type, "DRINKING_HOME");

homeEvnts = type.*events.times;
homeEvnts = homeEvnts(any(homeEvnts,2),:);
homeEvnts = homeEvnts(:,1); % times of stimulus onset

randEvnts = ~type.*events.times;
randEvnts = randEvnts(any(randEvnts,2),:);
randEvnts = randEvnts(:,1); % times of stimulus onset
%tf = isbetween(t,tlower,tupper);
end

