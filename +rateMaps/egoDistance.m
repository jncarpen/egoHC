function egoDistance(pos,SpikeTimes, refLoc, refLoc2, hd, sessNum)
%EGODISTANCE Summary of this function goes here
%   Detailed explanation goes here

% grab stuff
pos_ = pos{1,sessNum};
SpikeTimes_ = SpikeTimes;
hd_ = hd{1,sessNum};
hd_ = deg2rad(hd_)-pi; % range(-pi,+pi)


% Grab window limits for pos tracking
expansionFactor = 5;
xMin = nanmin(pos_(:,2))-expansionFactor;
xMax = nanmax(pos_(:,2))+expansionFactor;
yMin = nanmin(pos_(:,3))-expansionFactor;
yMax = nanmax(pos_(:,2))+expansionFactor;

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

% compute distance from rat to refLoc 
egoDist = sqrt((rlX-midX).^2 + (rlY-midY).^2);

% find time indices when cell spikes
idx = knnsearch(t, SpikeTimes_);
spkX = x1(idx); spkY = y2(idx);
spk_egoDist = egoDist(idx);

%% compute tuning curves
nBins = 20;
distEdges = linspace(0,max(egoDist),nBins);

%% egoDist
[spkEgoDistMap, mapAxis_EgoDist] = histcounts(spk_egoDist,distEdges);
[allEgoDistMap] = histcounts(egoDist,distEdges);

for i = 1:length(mapAxis_EgoDist)
    if i+1 <= length(mapAxis_EgoDist)
        binCtrs_egoDist(i) = ((mapAxis_EgoDist(i+1)-mapAxis_EgoDist(i))/2)+mapAxis_EgoDist(i);
    end
end
tcVals_egoDist = spkEgoDistMap./(allEgoDistMap*sampleRate + eps);

% smooth tuning curve values
tcVals_egoDist = imgaussfilt(tcVals_egoDist, 2, 'Padding', 'replicate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ref point #2
% break down 'refLoc2'
rlX2 = refLoc2(1,1);
rlY2 = refLoc2(1,2);

% compute distance from rat to refLoc 
egoDist2 = sqrt((rlX2-midX).^2 + (rlY2-midY).^2);

% find time indices when cell spikes
idx2 = knnsearch(t, SpikeTimes_);
spkX2 = x1(idx2); spkY = y2(idx2);
spk_egoDist2 = egoDist2(idx2);

%% compute tuning curves
nBins = 20;
distEdges2 = linspace(0,max(egoDist2),nBins);

%% egoDist
[spkEgoDistMap2, mapAxis_EgoDist2] = histcounts(spk_egoDist2,distEdges2);
[allEgoDistMap2] = histcounts(egoDist2,distEdges2);

for i = 1:length(mapAxis_EgoDist2)
    if i+1 <= length(mapAxis_EgoDist2)
        binCtrs_egoDist2(i) = ((mapAxis_EgoDist2(i+1)-mapAxis_EgoDist2(i))/2)+mapAxis_EgoDist2(i);
    end
end
tcVals_egoDist2 = spkEgoDistMap2./(allEgoDistMap2*sampleRate + eps);

% smooth tuning curve values
tcVals_egoDist2 = imgaussfilt(tcVals_egoDist2, 2, 'Padding', 'replicate');



%% plot
maximum = max(egoDist,egoDist2);
plot(binCtrs_egoDist, tcVals_egoDist, 'Color', 'k', 'LineWidth', 1.10)
hold on
plot(binCtrs_egoDist2, tcVals_egoDist2, 'LineStyle', ':','Color', 'r', 'LineWidth', 1.10)
title("egoDist")
xlabel("distance (cm)")
ylabel("fr (Hz)")
% xlim([0 maximum])
box off
end

