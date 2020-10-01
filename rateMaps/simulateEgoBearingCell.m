% SIMULATE EGOCENTRIC BEARING CELL:
% this script will use the data from Sebastian's OVC cell to
% simulate an egocentric bearing cell with (1) no noise and (2) with a
% uniform distribution of spontaneous firing noise.

% correct position data (cm) 
boxSize = 80; % in cm (according to Seb)
xLen = nanmax(positions(:,2))-nanmin(positions(:,2));
xLen2 = nanmax(positions(:,4))-nanmin(positions(:,4));
yLen = nanmax(positions(:,3))-nanmin(positions(:,3));
yLen2 = nanmax(positions(:,5))-nanmin(positions(:,5));
conFac_x = boxSize/xLen; conFac_y = boxSize/yLen;
conFac_x2 = boxSize/xLen2; conFac_y2 = boxSize/yLen2;

% parse corrected position data
t = positions(:,1); sampleRate = mode(diff(t));
x = (positions(:,2)-nanmin(positions(:,2)))*conFac_x; y = (positions(:,3)-nanmin(positions(:,3)))*conFac_y;
x2 = (positions(:,4)-nanmin(positions(:,4)))*conFac_x2; y2 = (positions(:,5)-nanmin(positions(:,5)))*conFac_y2;
posCorr = [t, x, y, x2, y2];

% calculate speed
x_smooth=medfilt1(x);y_smooth=medfilt1(y);
x2_smooth=medfilt1(x2);y2_smooth=medfilt1(y2);
speed_OVC = zeros(length(t), 2);
for i = 2:numel(x_smooth)-1
    speed_OVC(i, 1) = sqrt((x_smooth(i+1) - x_smooth(i-1))^2 + (y_smooth(i+1) - y_smooth(i-1))^2) / (t(i+1) - t(i-1));
    speed_OVC(i, 2) = sqrt((x2_smooth(i+1) - x2_smooth(i-1))^2 + (y2_smooth(i+1) - y2_smooth(i-1))^2) / (t(i+1) - t(i-1));
end
% pad the speed vectors
speed_OVC(1,1) = speed_OVC(2,1); speed_OVC(1,2) = speed_OVC(2,2);
speed_OVC(end, 2) = speed_OVC(end-1, 2); speed_OVC(end, 1) = speed_OVC(end-1, 1);


% find egocentric bearing for each timepoint
% rlX=42; rlY=37; % object
rlX=42; rlY=59; % object moved
midX=(x+x2)/2; midY=(y+y2)/2;
hd_OVC = rem(atan2d(y2-y, x2-x) + 180, 360);
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
% alloAng = rem(atan2d(midY-rlY, midX-rlX)+180, 360);
egoAng = alloAng - hd_OVC;

% define a range of egoAngles you are interested in,
% here I will use +/- 15 deg from 270 degrees. This means
% that the reference point is on the animal's RIGHT side.
angle_of_interest = 90; plus_minus = 15;
min_angle = angle_of_interest - plus_minus;
max_angle = angle_of_interest + plus_minus;

% find indices where egoAngle is within range
logical = egoAng>min_angle & egoAng<max_angle;
% logical = alloAng>min_angle & alloAng<max_angle;
idx = find(logical==1);

% select some of these indices to become spiketimes
throw_away = .5; % percentage to remove
sz = floor(length(idx)-length(idx)*throw_away); % size to keep
randIdx = datasample(idx, sz, 'Replace', false);
simTS = sort(t(randIdx), 'ascend'); % simulated timestamps

% speed threshold generated spikes
% speed threshold spike times
[STCorr] = speedThreshold_spikeTimes(posCorr, speed_OVC, simTS);


% find spike angle 
spk_index = knnsearch(t, STCorr);
spkAng = hd_OVC(spk_index);

%% plot
fig = figure('Position', [100 100 600 300]);
set(gcf,'color','w');
sgtitle("Simulated EBC")

% pathplot
subplot(1,3,1)
pathPlot_HD(posCorr, STCorr, hd_OVC);
xlim([nanmin(x), nanmax(x)])
ylim([nanmin(y), nanmax(y)])
pbaspect([1 1 1])
box off

% HD tuning curve
subplot(1,3,2)
tc_HD = analyses.turningCurve(spkAng, posCorr, sampleRate, 'smooth', 3, 'binWidth', 3);
% tc_HD = imgaussfilt(tc_HD, 2, 'Padding', 'circular'); % smooth
plot(tc_HD(:,1), tc_HD(:,2), 'Color', 'k', 'LineWidth', 1.10)
xlim([0 360])
title("HD tuning curve")
xlabel("head angle (deg)")
ylabel("fr (Hz)")
box off

% TC stats heatmap
% subplot(1,3,3)
% imagesc_output = tc_stats_heatmap(posCorr, hd_OVC, 0, STCorr, "False", "False", 'peak_from_shuff');
% imagesc(imagesc_output);
% set(gca,'YDir','normal')
% colormap(jet) % diverging
% colorbar
% % caxis([-1, 1])
% pbaspect([1 1 1])
% title('peak(real)-mean(peak(shuff))')
