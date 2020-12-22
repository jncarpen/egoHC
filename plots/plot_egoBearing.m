function [tcVals_egoAng, binCtrs_egoAng] = plot_egoBearing(position, ST, refLoc, doPlot)
%PLOT_EGOBEARING Thresholded!
% October 22, 2020
% This function can probably replace the other one (egoBearing.m)

% calculate head direction
[head_direction] = get_hd(position); % in deg

% Grab window limits for pos tracking
expansionFactor = 5;
xMin = nanmin(position(:,2))-expansionFactor;
xMax = nanmax(position(:,2))+expansionFactor;
yMin = nanmin(position(:,3))-expansionFactor;
yMax = nanmax(position(:,2))+expansionFactor;

% find point between 2 LEDs
t = position(:,1); 
sampleRate = mode(diff(t));
x1 = position(:,2); y1 = position(:,3);
x2 = position(:,4); y2 = position(:,5);
midX = (x1+x2)/2;
midY = (y1+y2)/2;

% break down 'refLoc'
rlX = refLoc(1,1);
rlY = refLoc(1,2);

% compute egocentric bearing
alloAng = atan2d(rlY-midY, rlX-midX)+180;
egoAng = alloAng - head_direction;
neg_idx = find(egoAng<0);
egoAng(neg_idx) = egoAng(neg_idx)+360;

%% apply speed threshold
threshold = 5; % cm/s
[P_logical, ST_thresh] = speed_threshold(position, ST, threshold);

% find time indices when cell spikes
idx = knnsearch(t, ST_thresh);
spkX = x1(idx); spkY = y2(idx);
spk_egoAng = egoAng(idx);

%% compute tuning curve
nBins = 40; % 9 degree bins
angEdges = linspace(0,360,nBins);
[spkEgoAngMap, mapAxis_EgoAng] = histcounts(spk_egoAng,angEdges);

% make the allEgoAngMap
for bin = 1:length(angEdges)-1
    in_this_bin = find(egoAng>angEdges(bin) & egoAng<=angEdges(bin+1));
    allEgoAngMap(1,bin) = sum(P_logical(in_this_bin));
end

for i = 1:length(mapAxis_EgoAng)
    if i+1 <= length(mapAxis_EgoAng)
        binCtrs_egoAng(i) = ((mapAxis_EgoAng(i+1)-mapAxis_EgoAng(i))/2)+mapAxis_EgoAng(i);
    end
end
tcVals_egoAng = spkEgoAngMap./(allEgoAngMap*sampleRate + eps); 

% smooth tuning curve values
tcVals_egoAng = imgaussfilt(tcVals_egoAng, 2, 'Padding', 'circular');

%% plot

switch doPlot
    case "True"
%         figure
        hold on;
        plot(binCtrs_egoAng, tcVals_egoAng, 'Color', [0 .4 .9], 'LineWidth', 1.10)
%         title("Egocentric Angle", "FontName", "Calibri Light", "FontSize", 20)
        ylabel("fr (Hz)"); xlabel("angle (deg)");
        xticks([90 180 270 360])
        box off
    case "False"
        % do nothing
end


end

