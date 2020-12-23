function [sim] = simulate_egoDist_cell(param)
%SIMULATE_EGO_CELL Simulate an egocentric bearing cell.
%   Inputs:
%   'pos_in'                    real position data [t x1 y1 x2 y2]
%   'ref_point'                 reference point [rlX, rlY]
%   'angle_of_interest'         range from 0-360 (need to wrap)
%   'plus_minus'                keep spikes within this range of
%                               angle_of_interest
%
%   Outputs:
%   'hd_sim'                    head direction values for simulated cell
%   'SpikeTimes_sim'            spiketimes for simulated cell
%   'SpikeTrain_sim'            spike train for simulated cell
%
% J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set box size- make this an input?
% boxSize = 80; (boxsize for seb's data)
% boxSize = 150; (boxsize for jan's data)

P = param.position;
ref_point = param.ref_point;
theta = param.theta;
radius = param.radius;

% parse position vector
t = P(:,1);
x = P(:,2); y = P(:,3);
x2 = P(:,4); y2 = P(:,5);

% get head_direction values
hd_sim = rem(atan2d(y2-y, x2-x) + 180, 360);

% parse ref_point
rlX = ref_point(1,1);
rlY = ref_point(1,2);

% define midpoint between two LEDs
midX=(x+x2)/2; midY=(y+y2)/2;

% calculate allocentric bearing
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);

% calculate egocentric bearing
egoAng = alloAng - hd_sim;

% correct for negative angles (egoAng)
neg_idx = find(egoAng<0);
egoAng(neg_idx) = egoAng(neg_idx)+360;

% define angles of interest
plus_minus_orien = 15; % how many degrees are acceptable
min_angle = theta - plus_minus_orien;
max_angle = theta + plus_minus_orien;

% calculate distance from reference point (for each timestamp)
dist_from_ref = sqrt((rlX - midX).^2+ (rlY - midY).^2);

% define distance range
plus_minus_distance = 15; % (in cm)
min_dist = radius - plus_minus_distance;
max_dist = radius + plus_minus_distance;

% find indices where egoAngle is within range
logical = egoAng>min_angle & egoAng<max_angle & dist_from_ref>min_dist & dist_from_ref<max_dist;
idx = find(logical==1); %logical = alloAng>min_angle & alloAng<max_angle;

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
figure
hold on;
set(gcf,'color','w');
pathPlot_hd(P, SpikeTimes_sim, hd_sim)
h1 = plot(rlX, rlY, 'o', 'MarkerSize', 8);
set(h1, 'markerfacecolor', 'k');
title("Egocentric Bearing + Distance")
% legend("path", "spikes", "refLoc", "Location", "northwestoutside")
hold off;

%% save things in a struct
sim.spiketrain = SpikeTrain_sim;
sim.spiketimes = SpikeTimes_sim;
sim.position = param.position;
sim.hd = hd_sim;
end