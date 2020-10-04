% CHECK_OVC:
% this script is just used for quickly visualizing the object-vector
% cell that Sebastian shared with me for the proof-of-concept analysis
% that we want to do on Jan Sigurd's data. J. Carpenter 2020.

% correct position data (cm) 
boxSize = 80; % in cm (according to Seb)
xLen = nanmax(positions(:,2))-nanmin(positions(:,2));
xLen2 = nanmax(positions(:,4))-nanmin(positions(:,4));
yLen = nanmax(positions(:,3))-nanmin(positions(:,3));
yLen2 = nanmax(positions(:,5))-nanmin(positions(:,5));
conFac_x = boxSize/xLen; conFac_y = boxSize/yLen;
conFac_x2 = boxSize/xLen2; conFac_y2 = boxSize/yLen2;

% parse corrected position data
t = positions(:,1);
x = (positions(:,2)-nanmin(positions(:,2)))*conFac_x; y = (positions(:,3)-nanmin(positions(:,3)))*conFac_y;
x2 = (positions(:,4)-nanmin(positions(:,4)))*conFac_x2; y2 = (positions(:,5)-nanmin(positions(:,5)))*conFac_y2;
posCorr = [t, x, y, x2, y2];
sampleRate = mean(diff(t));
STNow = cellTS; 

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

% speed threshold spike times
[STCorr] = speedThreshold_spikeTimes(posCorr, speed_OVC, STNow);


% calculate head direction (deg)
hd_OVC = rem(atan2d(y2-y, x2-x) + 180, 360);

% plot stuff
fig = figure('Position', [100 100 600 300]);
set(gcf,'color','w');
sgtitle("Object Trial")

% path plot (colored by HD)
subplot(2,2,1)
[imagesc_output, refVec] = tc_stats_heatmap(posCorr, hd_OVC, 0, STCorr, "False", "False", 'MVL_from_shuff');
pathPlot_HD(posCorr, STCorr, hd_OVC);
% hold on
for ii = 1:length(refVec)
plot(refVec(ii,1), refVec(ii,2), '.', 'Color', 'k')
hold on
end
yline(nanmin(y)); yline(nanmax(y));
xline(nanmin(x)); xline(nanmax(x));
xlim([nanmin(x), nanmax(x)])
ylim([nanmin(y), nanmax(y)])
pbaspect([1 1 1])
box off

% OVC Plot
subplot(2,2,2)
refCoord = [42,37]; % for object trial (cm)
% refCoord = [42, 59];
makeVecMaps(posCorr,STCorr,refCoord);
box off

% HD tuning curve
subplot(2,2,3)
spkang = getSpikeAngle(hd_OVC, STCorr);
tc = analyses.turningCurve(spkang, hd_OVC, sampleRate,'smooth', 2.5, 'binWidth', 5);
plot(tc(:,1), tc(:,2), 'Color', 'k', 'LineWidth', 1.10)
xlim([0 360])
xlabel("head angle (deg)")
ylabel("fr (hz)")
box off

% egoBearing TC stats
subplot(2,2,4)
imagesc(imagesc_output);
set(gca,'YDir','normal')
colormap(jet) % diverging
colorbar
caxis([-1, 1])
pbaspect([1 1 1])
title('MVL(real)-mean(MVL(shuff))')

