function [speed_thresh_spikeTimes] = speedThreshold_spikeTimes(pos_, speed_cm, SpikeTimes_)
%SPEEDTHRESHOLD_SPIKETRAIN Summary of this function goes here
%   Detailed explanation goes here

% set conversion factor (units to cm)
% I need to go back and correct the position data in matlab, but for
% now this will do

% grab speed for first LED
spdLED1 = speed_cm(:,1);
minSpd = nanmin(spdLED1);
maxSpd = nanmax(spdLED1);
t = pos_(:,1);
sampleRate = mode(diff(t));

% find closest timestamp & store values
idx = knnsearch(t, SpikeTimes_);
spkSpd = spdLED1(idx);
greater_than_5_idx = find(spkSpd>5);
speed_thresh_spikeTimes = SpikeTimes_(greater_than_5_idx);
end

