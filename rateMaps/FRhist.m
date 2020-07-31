function [out] = FRhist(pos, SpikeTrain, speed)
%FRHIST Bin SpikeTrain into 200 ms bins and make histogram.
%
%   With a sampling frequency of 0.020 seconds, we want to combine every 10
%   bins to yield 200 ms bins.
%
%   USAGE
%   [out] = FRhist(SpikeTimes, pos, <options>);
%
%   INPUT
%   SpikeTrain:         Binned SpikeTimes
%   pos:                [t x y x2 y2]
%
%
%   OUTPUT
%   out:                lorem
%
%
%   =========================================================================
%     Properties    Values
%   -------------------------------------------------------------------------
%     'speedFilt'   determines whether a firing rate histogram is created for
%                   all spikes (0), for spikes that fire when the speed of the
%                   animal is LESS than 5cm/s (1) or when speed>5cm/s (2).
%                   Default is to create a histogram for *all* cells
%                   ('speedFilt' = 0).
%   =========================================================================
%
%   Jordan Carpenter, 2020.

%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%

% parse inputs
% p = inputParser;
% addRequired(p,'SpikeTimes',@isnumeric)
% addRequired(p,'pos',@isnumeric)
% addParameter(p,'speedFilt',[0],@isinteger)% name-value pair
% addOptional(p, 'speed', @isnumeric)
% parse(p,varargin{:});

% determine whether a speed threshold should be applied
% if speedFilt == 1
%     
% elseif speedFilt == 2
% end


% pull out important stuff
t = pos(:,1); 
minTime = nanmin(t);
maxTime = nanmax(t);
binWidth = 0.2; % 200 ms bins

% build histogram
edges = minTime:binWidth:maxTime; 
[N] = histcounts(SpikeTimes,edges); % N = number of occurances/bin

% compute average FR for each bin (Hz)
firingRate = N./binWidth;

% plot histogram of average firing rates
nBins = 100; % # of bins
FRhist = histcounts(firingRate, 100);


end

