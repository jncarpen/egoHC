function [sim] = simulate_place_depreciated(param)
%SIMULATE_PLACE:
%   Simulate a place cell with given parameters.
%   INPUT STRUCT-
%   param.position:
%   param.ctr_mass:
%   param.r:
%
%   OUTPUT STRUCT-
%   sim.spiketimes:
%   sim.spiketrain:
%   sim.hd:
%
% J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse position vector
P = param.position;
t = P(:,1);
x = P(:,2); y = P(:,3);
x2 = P(:,4); y2 = P(:,5);

% get head_direction values
head_direction = atan2d(y2-y, x2-x) + 180;

% parse center of mass
ctr_massX = param.ctr_mass(1,1);
ctr_massY = param.ctr_mass(1,2);

% define midpoint between two LEDs
midX=(x+x2)/2; midY=(y+y2)/2;

% calculate distance from place field center [of mass]
dist_from_ctr = sqrt((ctr_massX - midX).^2+ (ctr_massY - midY).^2);

% define distance range
max_distance = 5; % (in cm)
min_dist = 0;
max_dist = param.r + max_distance;

% find indices where egoAngle is within range
logical = dist_from_ctr>min_dist & dist_from_ctr<max_dist;
idx = find(logical==1); %logical = alloAng>min_angle & alloAng<max_angle;

% define how many spikes to keep
throw_away = .45;
sz = floor(length(idx)-length(idx)*throw_away);
randIdx = datasample(idx, sz, 'Replace', false);
foreground_spikes = t(randIdx);


% number of background spikes to add (noise)
len_background = sz*.3; 
% need to add the noise in this line

% sort the simulated spikes
simTS = sort(foreground_spikes, 'ascend'); % simulated timestamps

% speed threshold spikes (keep >5cm/s)
x_smooth=medfilt1(x);y_smooth=medfilt1(y);
speed_OVC = zeros(length(t), 1);
for i = 2:numel(x_smooth)-1
    speed_OVC(i) = sqrt((x_smooth(i+1) - x_smooth(i-1))^2 + (y_smooth(i+1) - y_smooth(i-1))^2) / (t(i+1) - t(i-1));
end

% make a spike train
stopTime = nanmax(t); startTime = nanmin(t);
SpikeTimes_sim = simTS(simTS < stopTime & simTS > startTime); % remove times outside of recording
edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate
binnedSpikes = histcounts(SpikeTimes_sim,edgesT); % bin those spikes baby
speed_idx = find(speed_OVC<5); % find indices when animal was moving slow
binnedSpikes(speed_idx)=0; % get rid of spikes when animal was moving slow
binnedSpikes = imgaussfilt(binnedSpikes, 2, 'Padding', 'replicate'); % smooth ST
SpikeTrain_sim = binnedSpikes;

% find position stamps for the spike times
spkIdx = knnsearch(t, SpikeTimes_sim);
xKeep = x(spkIdx); yKeep = y(spkIdx);
[N,Xedges,Yedges] = histcounts2(xKeep,yKeep);




% show user their cell!
figure
hold on;
set(gcf,'color','w');
pathPlot_hd_Gaussian(P, SpikeTimes_sim, head_direction);
title("Place Cell")
hold off;

% save information in a struct
sim.spiketimes = SpikeTimes_sim;
sim.spiketrain = SpikeTrain_sim;
sim.hd = head_direction;


end

