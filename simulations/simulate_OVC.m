function [SpikeTimes_sim, SpikeTrain_sim, hd_sim] = simulate_OVC(pos_in, ref_point, angle_of_interest, radius)
%SIMULATE_OVC Summary of this function goes here
%   Detailed explanation goes here
%   Inputs:
%   'pos_in'                    real position data [t x1 y1 x2 y2]
%   'ref_point'                 reference point [rlX, rlY]
%   'angle_of_interest'         range from 0-360 (need to wrap)
%   'plus_minus'                keep spikes within this range of
%                               angle_of_interest
%
%   Outputs:
%   'pos_sim'                   position data for simulated cell
%   'hd_sim'                    head direction values for simulated cell
%   'SpikeTimes_sim'            spiketimes for simulated cell
%   'SpikeTrain_sim'            spike train for simulated cell
%
% J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse position vector
t = pos_in(:,1);
x = pos_in(:,2); y = pos_in(:,3);
x2 = pos_in(:,4); y2 = pos_in(:,5);

% get head_direction values
hd_sim = rem(atan2d(y2-y, x2-x) + 180, 360);

% parse ref_point
rlX = ref_point(1,1);
rlY = ref_point(1,2);

% define midpoint between two LEDs
midX=(x+x2)/2; midY=(y+y2)/2;

% calculate allocentric bearing
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);

% define angles of interest
plus_minus_orien = 10; % how many degrees are acceptable
min_angle = angle_of_interest - plus_minus_orien;
max_angle = angle_of_interest + plus_minus_orien;

% calculate distance from reference point (for each timestamp)
dist_from_ref = sqrt((rlX - midX).^2+ (rlY - midY).^2);

% define distance range
plus_minus_distance = 5; % (in cm)
min_dist = radius - plus_minus_distance;
max_dist = radius + plus_minus_distance;

% find indices where egoAngle is within range
logical = alloAng>min_angle & alloAng<max_angle & dist_from_ref>min_dist & dist_from_ref<max_dist;
idx = find(logical==1); %logical = alloAng>min_angle & alloAng<max_angle;

% define how many spikes to keep
throw_away = .35;
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
h1 = plot(rlX, rlY, 'o', 'MarkerSize', 8);
set(h1, 'markerfacecolor', 'k');
% legend("path", "spikes", "refLoc", "Location", "northwestoutside")
title("Allocentric Bearing + Distance (OVC)")
hold off;


end

