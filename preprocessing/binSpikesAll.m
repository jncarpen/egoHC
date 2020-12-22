function spikeTrain = binSpikesAll(trackingtimes, cellSpikes)
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

% define sample rate- use for speed stuff later ** 
% sampleRate = mode(diff(trackingtimes));

% intiialize empty cell array to store smoothed spike info for each cell
spikeTrain = cell(1, length(cellSpikes));


% loop through all cells
for neuron = 1:length(cellSpikes)
    
    % get spikes for current cell
    S = cellSpikes{1,neuron};
    
    % remove spike times that are outside the range of tracking times
    S = S(S < stopTime & S > startTime);
    
    % bin spike data (by time)
    edgesT = linspace(startTime,stopTime,numel(trackingtimes)+1); % binsize is close to video frame rate
    
    % histcount of spikes per time bin
    binnedSpikes = histcounts(S,edgesT);
    
    % smooth spike train
    sigma = 2;
    binnedSpikes = imgaussfilt(binnedSpikes, sigma, 'Padding', 'replicate');
    
    % save to cell array
    spikeTrain{1,neuron} = binnedSpikes;
    
    
%% Smooth spikes:
%     % initialize bin size for gaussian window
%     binSz = 5;
%     
%     % smooth time-binned data w/ a gaussian filter
%     binnedSpikes_smooth = smoothdata(binnedSpikes, 'gaussian', binSz);
%     
%     % store smoothed values for each neuron in a cell array
%     spikeTrain{1, neuron} = binnedSpikes_smooth;

end

