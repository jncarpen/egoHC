function [spikeTrain, spikeTrain_smooth] = binSpikes(trackingtimes, S)
%BINSPIKES Summary of this function goes here
%   Input cell array of cellSpikes for a single session

%   Output is the *smoothed spike train (can adjust the width of the
%   gaussian window)

%   Everything should be in SECONDS **

% if trackingtimes is in ms;convert to seconds
% trackingtimes = trackingtimes/1000;

% compute first and last timestamp
startTime = trackingtimes(1);
stopTime = trackingtimes(end);

% remove spike times that are outside the range of tracking times
S = S(S < stopTime & S > startTime);

% bin spike data (by time)
edgesT = linspace(startTime,stopTime,numel(trackingtimes)+1);

% histcount of spikes per time bin
spikeTrain = histcounts(S, edgesT);

% smooth spiketrain
sigma = 2;
spikeTrain_smooth = imgaussfilt(spikeTrain, sigma, 'Padding', 'replicate');

end

