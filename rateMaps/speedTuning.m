function [] = speedTuning(speed)
%SPEEDTUNING Compute speed tuning curve
%   This function isn't done yet

maxSpeed = max(speed);
minSpeed = min(speed);
numBins = %?;

speedBins = linspace(minSpeed,maxSpeed,numBins+1);
speedTC = nan(numBins,1); % speed tuningcurve _init_
end

