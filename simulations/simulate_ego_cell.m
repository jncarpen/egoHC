function [sim] = simulate_ego_cell(param)
%SIMULATE_EGO_CELL Simulate an egocentric bearing cell.
%   INPUT STRUCT -
%   param.position
%   param.ref_point
%   param.theta
%   param.hd
%
%   OUTPUT STRUCT - 
%   sim.hd                      head direction values for simulated cell
%   sim.spiketimes              spiketimes for simulated cell
%   sim.spiketrain              spike train for simulated cell
%
% J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse position vector
t = param.P(:,1);
x = param.P(:,2); y = param.P(:,3);
hd_sim = param.hd;

% parse ref_point
rlX = param.ref_point(1,1);
rlY = param.ref_point(1,2);

% define midpoint between two LEDs
midX = x; midY = y;

% calculate allocentric bearing
egoAng = mod((atan2d(rlY-midY, rlX-midX)+180) - hd_sim, 360);

% define angles of interest
plus_minus_orien = 30; % how many degrees are acceptable
min_angle = param.theta - plus_minus_orien;
max_angle = param.theta + plus_minus_orien;

% find indices where egoAngle is within range
logical = egoAng>min_angle & egoAng<max_angle;
idx = find(logical==1);

% define how many spikes to keep
throw_away = .01;
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

sim.ST = SpikeTimes_sim;
sim.spiketrain = SpikeTrain_sim;

% show user their simulated cell!
% figure
% hold on;
% set(gcf,'color','w');
% pathPlot_hd(param.position, SpikeTimes_sim, hd_sim) % from -180 to +180
% c1 = colorbar; c1.Ticks = [0 90 180 270 360]; c1.FontSize = 25;
% h1 = plot(rlX, rlY, 'o', 'MarkerSize',15);
% set(h1, 'markerfacecolor', 'k');
% color_plot_title = strcat('theta_{pref} = ', sprintf('%.f',param.theta));
% title(color_plot_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
% 
% 
% figure
% pathPlot_quiver(param.position, SpikeTimes_sim, hd_sim)
% 
% 
% figure
% map = analyses.map(param.position, SpikeTimes_sim, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
% peakRate = nanmax(nanmax(map.z));
% rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
% plot.colorMap(map.z)
% pbaspect([1 1 1])
% colormap(gca,'jet')
% c2 = colorbar; c2.FontSize = 25;
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
% box off

end
