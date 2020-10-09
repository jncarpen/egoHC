function [allSpdCnts, hiSpdCnts, loSpdCnts, histEdges] = FRhist(pos, SpikeTimes, speed)
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

% pull out important stuff
t = pos(:,1); 
minTime = nanmin(t);
maxTime = nanmax(t);
binWidth = 0.2; % 200 ms bins

% get spike speed
spdLED1 = speed(:,1);
idx = knnsearch(t, SpikeTimes);
spkSpd = spdLED1(idx);

% get indices for times where speed was < or > 5 cm/s
over5idx = find(spkSpd>5);
under5idx = find(spkSpd<5);

ST_over5 = SpikeTimes(over5idx);
ST_under5 = SpikeTimes(under5idx);

%% ALL FIRING RATES

% build histogram
edges = minTime:binWidth:maxTime; 
[N] = histcounts(SpikeTimes,edges); % N = number of occurances/bin

% compute average FR for each bin (Hz)
firingRateALL = N./binWidth;

% plot histogram of average firing rates
% edgesFR = linspace(0.001,50,50);
edgesFR = logspace(-3, 1.7,50); % fom 0.001 to ~50
%edgesFR = logspace(-3,1.30103, 31); % from 0.0010 to 20.0000
[allSpdCnts, histEdges] = histcounts(firingRateALL, edgesFR);

%% OVER 5 CM/S

% build histogram
edges = minTime:binWidth:maxTime; 
[N_over5] = histcounts(ST_over5,edges); % N = number of occurances/bin

% compute average FR for each bin (Hz)
firingRateHI = N_over5./binWidth;

% plot histogram of average firing rates
[hiSpdCnts, ~] = histcounts(firingRateHI, edgesFR);

%% UNDER 5 CM/S

% build histogram
edges = minTime:binWidth:maxTime; 
[N_under5] = histcounts(ST_under5,edges); % N = number of occurances/bin

% compute average FR for each bin (Hz)
firingRateLO = N_under5./binWidth;

% plot histogram of average firing rates
[loSpdCnts, ~] = histcounts(firingRateLO, edgesFR);

end

