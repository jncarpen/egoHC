function [SpikeTimes_sim, SpikeTrain_sim, hd_sim] = simulate_HD(pos_in, angle_of_interest)
%SIMULATE_HD Simulate a head-direction cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% parse position vector
t = pos_in(:,1);
x = pos_in(:,2); y = pos_in(:,3);
x2 = pos_in(:,4); y2 = pos_in(:,5);

% get head_direction values
hd_sim = rem(atan2d(y2-y, x2-x) + 180, 360);


% define angles of interest
plus_minus_orien = 5; % +/- orientation
min_angle = angle_of_interest - plus_minus_orien;
max_angle = angle_of_interest + plus_minus_orien;

% find indices where egoAngle is within range
logical = hd_sim>min_angle & hd_sim<max_angle;
idx = find(logical==1);

% define how many spikes to keep
throw_away = .25;
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

% show user their cell!
% figure
hold on;
set(gcf,'color','w');
pathPlot_hd(pos_in, SpikeTimes_sim, hd_sim)
title("Head-Direction")
hold off;

end

