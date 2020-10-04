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
xCorr = posCorr(:,2); yCorr = posCorr(:,3);

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
angle_of_interest = 60; plus_minus = 20;
min_angle = angle_of_interest - plus_minus;
max_angle = angle_of_interest + plus_minus;

% find indices where egoAngle is within range
logical = egoAng>min_angle & egoAng<max_angle;
% logical = alloAng>min_angle & alloAng<max_angle;
idx = find(logical==1);

% select some of these indices to become spiketimes
throw_away = .2; % percentage to remove
sz = floor(length(idx)-length(idx)*throw_away); % size to keep
len_background = sz*.3; % number of background spikes to add

% sample (1-throw_away)*100 (percent) of the possible values
randIdx = datasample(idx, sz, 'Replace', false);
foreground_spikes = t(randIdx);

% generate a uniform distribution of background spikes
% background_spikes = (nanmax(t)-nanmin(t)).*rand(floor(len_background),1) + nanmin(t);

% concatenate fore/background and sort in ascending order
% simTS = sort([foreground_spikes;background_spikes], 'ascend'); % simulated timestamps
simTS = sort(foreground_spikes, 'ascend'); % simulated timestamps

% speed threshold generated spikes
[STCorr] = speedThreshold_spikeTimes(posCorr, speed_OVC, simTS);

% find spike angles
spk_index = knnsearch(t, STCorr);
spkAng = hd_OVC(spk_index);
spkBearing = egoAng(spk_index);

% find spike positions
xSpk=xCorr(spk_index); ySpk=yCorr(spk_index);
r = 1; % circle radius
horiz = r*cos(spkBearing); vert = r*sin(spkBearing);

%% Quiver pathplot
% identify data in *each quadrant*
figure(2)
title("Egocentric Bearing (60 deg)")
pbaspect([1 1 1])
xlim([nanmin(x), nanmax(x)])
ylim([nanmin(y), nanmax(y)])
box off

% choose a theta
theta = spkBearing;

% Data is organized as (x, y, theta in degrees)
data = [xSpk, ySpk, theta];

% Identify data from 1st quadrant and store in q1, identify data from 2nd quadrant
% goes to q2, and so on. A value of 1 indicates it is in the quadrant, 0 otherwise
q1 = (data(:,3) <= 90);
q2 = (data(:,3) > 90) .* (data(:,3) <= 180);
q3 = (data(:,3) > 180) .* (data(:,3) <= 270);
q4 = (data(:,3) > 271) .* (data(:,3) <= 360);
hold on;

% Use QUIVER to specify the start point (tail of the arrow) and direction based on angle
% q1, q2, q3, and q4 are used to generate four different QUIVER handles (h1, h2, h3, and h4)
% This is necessary for varying colors based on direction
% Based on equations: x = x0 + r*cos(theta), y = y0 + r*sin(theta)
% In the usage below, x0 = data(:,1), y0 = data(:,2), theta = data(:,3) * pi / 180
% Can also specify a scale factor as the last argument to quiver (not specified below)
h1 = quiver(data(q1 == 1,1), data(q1 == 1,2), cos(data(q1 == 1,3) * pi/180), sin(data(q1 == 1,3) * pi/180));
h2 = quiver(data(q2 == 1,1), data(q2 == 1,2), cos(data(q2 == 1,3) * pi/180), sin(data(q2 == 1,3) * pi/180)); % sin is negative in 2nd quadrant
h3 = quiver(data(q3 == 1,1), data(q3 == 1,2), cos(data(q3 == 1,3) * pi/180), sin(data(q3 == 1,3) * pi/180));
h4 = quiver(data(q4 == 1,1), data(q4 == 1,2), cos(data(q4 == 1,3) * pi/180), sin(data(q4 == 1,3) * pi/180)); % cos is negative in 4th quadrant

% Set colors to red for 1st quadrant, blue for 2nd, green for 3rd, cyan for 4th
% Also, turn scaling off. get(h1) will return additional property-value pairs
set(h1, 'Color', 'k', 'AutoScale', 'off')


% set(h1, 'Color', [.8 0 .5], 'AutoScale', 'off')
% set(h2, 'Color', [.9 .5 .3], 'AutoScale', 'off')
% set(h3, 'Color', [.1 0 .8], 'AutoScale', 'off')
% set(h4, 'Color', [0 .75 .5], 'AutoScale', 'off')
% Done plotting
hold off;



% figure
% set(gcf,'color','w');
% 
% subplot(1,2,1)
% pathPlot_HD(posCorr, STCorr, hd_OVC);
% xlim([nanmin(x), nanmax(x)])
% ylim([nanmin(y), nanmax(y)])
% pbaspect([1 1 1])
% box off
% 
% subplot(1,2,2)
% quiver(xSpk, ySpk, horiz, vert)
% xlim([nanmin(x), nanmax(x)])
% ylim([nanmin(y), nanmax(y)])
% pbaspect([1 1 1])
% box off

%% plot
fig = figure('Position', [100 100 600 300]);
sgtitle("Simulated EBC")

% pathplot
subplot(3,2,6)
pathPlot_HD(posCorr, STCorr, hd_OVC);
title("270 deg")
xlim([nanmin(x), nanmax(x)])
ylim([nanmin(y), nanmax(y)])
pbaspect([1 1 1])
box off

quiver(posCorr(:,2), posCorr(:,3), hdOVC)

% reference point plot
% subplot(2,2,2)
% [imagesc_output, refVec] = tc_stats_heatmap(posCorr, hd_OVC, 0, STCorr, "False", "False", 'MVL_from_shuff');
% for ii = 1:length(refVec)
% plot(refVec(ii,1), refVec(ii,2), '.')
% hold on
% end
% yline(nanmin(y)); yline(nanmax(y));
% xline(nanmin(x)); xline(nanmax(x));
% pbaspect([1 1 1])
% box off

% HD tuning curve
subplot(2,2,3)
tc_HD = analyses.turningCurve(spkAng, posCorr, sampleRate, 'smooth', 3, 'binWidth', 3);
% tc_HD = imgaussfilt(tc_HD, 2, 'Padding', 'circular'); % smooth
plot(tc_HD(:,1), tc_HD(:,2), 'Color', 'k', 'LineWidth', 1.10)
xlim([0 360])
title("HD tuning curve")
xlabel("head angle (deg)")
ylabel("fr (Hz)")
box off

% MVL/TC Stats Map
% subplot(2,2,4)
% imagesc(imagesc_output);
% set(gca,'YDir','normal')
% colormap(jet) % diverging
% colorbar
% caxis([-1, 1])
% pbaspect([1 1 1])
% title('MVL(real)-mean(MVL(shuff))')
