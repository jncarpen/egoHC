function [TC_taskPhz] = taskPhz_tuningCurve(pos, SpikeTrain, taskPhz)
%TASKPHZ_TUNINGCURVE Summary of this function goes here
%   Detailed explanation goes here
t = pos(:,1);
Fs = mode(diff(t));
bins = [0 1 2 3 4 5 6];

% Find TP at time of spikes
idx = find(SpikeTrain>=1);
spkPhz = taskPhz(idx);

[phzOcc,~] = histcounts(taskPhz, bins);

% spikes in each task phase
spkPerPhz = histcounts(spkPhz, bins);

% compute tuning curve values
TC_taskPhz = spkPerPhz./phzOcc * Fs;

% smooth tuning curve values
TC_taskPhz = imgaussfilt(TC_taskPhz, 2, 'Padding', 'replicate');

end

