function [spkhd] = getSpikeAngle(hd, SpikeTrain)
%UNTITLED 
%   USAGE

% find indices where cell spiked
idx = find(SpikeTrain >= 1);

% pull out hd value at that index (in DEG)
spkhd = hd(idx);

end

