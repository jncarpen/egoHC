function [sim] = simulate_egoDist_cell(param)
%SIMULATE_EGODIST_CELL
%   INPUTS-
%   param.theta = 0; % facing toward object
%   param.Z = angular variable;
%   param.P = P;
%   param.rp = [75, 75];
%   param.radius = 10;
%   J. Carpenter, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unpack root structure
theta = param.theta;
radius = param.radius;
Z = param.Z;
rlX = param.rp(1);
rlY = param.rp(2);

t = param.P(:,1);
[~,c] = size(param.P);
switch c
    case 3
        x1 = param.P(:,2); y1 = param.P(:,3);
        midX = x1; midY = y1;
    case 5
        x1 = param.P(:,2); y1 = param.P(:,3);
        x2 = param.P(:,4); y2 = param.P(:,5);
        midX=(x1+x2)/2; midY=(y1+y2)/2;
end

% calculate allocentric bearing
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);

% calculate egocentric bearing
egoAng = mod(alloAng - Z, 360);

% define angles of interest
plus_minus_orien = 35; % how many degrees are acceptable
min_angle = theta - plus_minus_orien;
max_angle = theta + plus_minus_orien;

% calculate distance from reference point (for each timestamp)
dist_from_ref = sqrt((rlX - midX).^2+ (rlY - midY).^2);

% define distance range
plus_minus_distance = 20; % (in cm)
min_dist = radius - plus_minus_distance;
max_dist = radius + plus_minus_distance;

% find indices where egoAngle is within range
logical = egoAng>min_angle & egoAng<max_angle & dist_from_ref>min_dist & dist_from_ref<max_dist;
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
x_smooth=medfilt1(x1);y_smooth=medfilt1(y1);
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
pathPlot_hd(param.P, SpikeTimes_sim, Z)
h1 = plot(rlX, rlY, 'o', 'MarkerSize', 8);
set(h1, 'markerfacecolor', 'k');
title("Egocentric Bearing + Distance")
% legend("path", "spikes", "refLoc", "Location", "northwestoutside")
hold off;

%% save things in a struct
sim.TR = SpikeTrain_sim;
sim.ST = SpikeTimes_sim;
end
