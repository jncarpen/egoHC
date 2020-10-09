function [spkPos,spkInd, posCounts] = spikePos(SpikeTimes, pos)
%SPIKEPOS Summary of this function goes here
%   Detailed explanation goes here

t = pos(:,1); % time
XY = pos; % duplicate pos matrix
XY(:,1) = []; % remove time column
posCounts = zeros(size(t));

% sampleRate = mode(diff(t));
% minTime = nanmin(t);
% maxTime = nanmin(t) + sampleRate;

% make sure SpikeTimes is a column vector
SpikeTimes = SpikeTimes(:);

N = size(SpikeTimes, 1);
spkPos = zeros(N, size(pos, 2)); % N rows, 5 columns
spkInd = zeros(N, 1);

% Find temporal spike index (TSI) by finding k-nearest neighbors
TSI = knnsearch(t, SpikeTimes); 

count = 0;
for i = length(TSI)
    count = count + 1;
    ind = TSI(i);
    spkPos(count, 1) = t(i);
    spkPos(count, 2:3) = XY(ind, 2:3);
    spkInd(count) = ind;
    posCounts(ind) = posCounts(ind) + 1;
end

% remove empty rows
% spkPos = spkPos(1:count, :); 
% spkInd = spkInd(1:count);

end

