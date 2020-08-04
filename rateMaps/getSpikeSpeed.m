function [spkSpd] = getSpikeSpeed(speed,SpikeTimes)
%GETSPIKESPEED Find speed of the animal each time the cell fired.
%   INPUT
%   speed:              [spdLED1 spdLED2]
%   SpikeTimes:         1xS vector, where S is number of times the cell fired.
%   
%   OUTPUT
%   spkSpd:             1xS vector (should be same length as SpikeTimes)
%                       with the speed value (in cm/s) of the animal at each timestamp that the
%                       cell fired.
%
%   Jordan Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab speed for first LED
spdLED1 = speed(:,1);

% find closest timestamp & store values
idx = knnsearch(spdLED1, SpikeTimes);
spkSpd = spdLED1(idx);

end

