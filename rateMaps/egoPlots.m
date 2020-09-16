function egoPlots(pos,SpikeTimes, refLoc, hd, sessInfo, sessNum, unitNum)
%EGOPLOTS Summary of this function goes here
%   'pos':         [t x1 y1 x2 y2]
%   'refLoc':      [x1 y1]
%
%   Notes:
%   Convert distance to centimeters and find out what the reference is for
%   head direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab stuff
pos_ = pos{1,sessNum};
SpikeTimes_ = SpikeTimes{1,sessNum}{1,unitNum};
hd_ = hd{1,sessNum};
hd_ = deg2rad(hd_)-pi; % range(-pi,+pi)

% Grab window limits for pos tracking
xMin = str2double(sessInfo{1,sessNum}.window_min_x{1,1})-10;
xMax = str2double(sessInfo{1,sessNum}.window_max_x{1,1})+10;
yMin = str2double(sessInfo{1,sessNum}.window_min_y{1,1})-10;
yMax = str2double(sessInfo{1,sessNum}.window_max_y{1,1})+10;

% find point between 2 LEDs
t = pos_(:,1); 
sampleRate = mode(diff(t));
x1 = pos_(:,2); y1 = pos_(:,3);
x2 = pos_(:,4); y2 = pos_(:,5);
midX = (x1+x2)/2;
midY = (y1+y2)/2;

% break down 'refLoc'
rlX = refLoc(1,1);
rlY = refLoc(1,2);

% compute egocentric angle to reference loc
% which one is correct?
egoAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
% egoAng = rem(atan2d(midY-rlY, midX-rlX)+180, 360);
egoAng = deg2rad(egoAng)-pi; % range(-pi,+pi)

% compute distance from rat to refLoc 
egoDist = sqrt((rlX-midX).^2 + (rlY-midY).^2);

% find time indices when cell spikes
idx = knnsearch(t, SpikeTimes_);
spkX = x1(idx); spkY = y2(idx);
spk_egoAng = egoAng(idx);
spk_egoDist = egoDist(idx);
spk_hd = hd_(idx);

%% compute tuning curves
nBins = 20;
angEdges = linspace(-pi,pi,nBins); 
distEdges = linspace(0,max(egoDist),nBins);

%% egoAng
[spkEgoAngMap, mapAxis_EgoAng] = histcounts(spk_egoAng,angEdges);
[allEgoAngMap] = histcounts(egoAng,angEdges);

for i = 1:length(mapAxis_EgoAng)
    if i+1 <= length(mapAxis_EgoAng)
        binCtrs_egoAng(i) = ((mapAxis_EgoAng(i+1)-mapAxis_EgoAng(i))/2)+mapAxis_EgoAng(i);
    end
end
tcVals_egoAng = spkEgoAngMap./(allEgoAngMap*sampleRate + eps); 

%% egoDist
[spkEgoDistMap, mapAxis_EgoDist] = histcounts(spk_egoDist,distEdges);
[allEgoDistMap] = histcounts(egoDist,distEdges);

for i = 1:length(mapAxis_EgoDist)
    if i+1 <= length(mapAxis_EgoDist)
        binCtrs_egoDist(i) = ((mapAxis_EgoDist(i+1)-mapAxis_EgoDist(i))/2)+mapAxis_EgoDist(i);
    end
end
tcVals_egoDist = spkEgoDistMap./(allEgoDistMap*sampleRate + eps);

%% hd
[spkHDMap, mapAxis_HD] = histcounts(spk_hd,angEdges);
[allHDMap] = histcounts(hd_,angEdges);

for i = 1:length(mapAxis_HD)
    if i+1 <= length(mapAxis_HD)
        binCtrs_HD(i) = ((mapAxis_HD(i+1)-mapAxis_HD(i))/2)+mapAxis_HD(i);
    end
end
tcVals_HD = spkHDMap./(allHDMap*sampleRate + eps);


%% plot
figure('Position', [100 100 700 300], 'Color','white');
subplot(2,3,1)
plot(x1, y1, 'Color', [.7 .7 .7])
hold on
scatter(spkX, spkY, [30], spk_egoAng, '.')
axis square
box off
title("egoAng")
colorbar(gca, "hsv", 'eastoutside')
set(gca,'xtick',[],'ytick',[])
xlim([xMin xMax])
ylim([yMin yMax])

subplot(2,3,2)
plot(x1, y1, 'Color', [.7 .7 .7])
hold on
scatter(spkX, spkY, [30], spk_egoDist, '.')
axis square
box off
title("egoDist")
colorbar(gca, "hsv", 'eastoutside')
set(gca,'xtick',[],'ytick',[])
xlim([xMin xMax])
ylim([yMin yMax])

subplot(2,3,3)
plot(x1, y1, 'Color', [.7 .7 .7])
hold on
scatter(spkX, spkY, [30], spk_hd, '.')
axis square
box off
title("HD")
colorbar(gca, "hsv", 'eastoutside')
set(gca,'xtick',[],'ytick',[])
xlim([xMin xMax])
ylim([yMin yMax])

subplot(2,3,4)
plot(binCtrs_egoAng, tcVals_egoAng, 'Color', 'k', 'LineWidth', 1.10)
title("egoAng")
xlabel("angle (rad)")
ylabel("fr (Hz)")
xlim([-pi pi])

subplot(2,3,5)
plot(binCtrs_egoDist, tcVals_egoDist, 'Color', 'k', 'LineWidth', 1.10)
title("egoDist")
xlabel("angle (rad)")
ylabel("fr (Hz)")
xlim([0 max(egoDist)])

subplot(2,3,6)
plot(binCtrs_HD, tcVals_HD, 'Color', 'k', 'LineWidth', 1.10)
title("HD")
xlabel("angle (rad)")
ylabel("fr (Hz)")
xlim([-pi pi])

return


end

