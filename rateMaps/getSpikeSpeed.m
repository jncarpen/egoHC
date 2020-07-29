function [spkSpd] = getSpikeSpeed(speed,SpikeTrain)
%GETSPIKESPEED Summary of this function goes here
%   Detailed explanation goes here

% find indices where cell spiked
idx = find(SpikeTrain >= 1);

% pull out hd value at that index (in DEG)
spkSpd = speed(idx);

end

