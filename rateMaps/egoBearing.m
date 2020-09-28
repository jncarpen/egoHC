function [tcVals_egoAng] = egoBearing(pos_cm, SpikeTimes, refLoc, refLoc2, hd, sessNum, doPlot, deg_or_rad)
%EGOBEARING Summary of this function goes here

% grab stuff
pos_ = pos_cm{1,sessNum};
SpikeTimes_ = SpikeTimes;
hd_ = hd{1,sessNum};

% decide whether to convert to radians or keep in degrees
if deg_or_rad == "rad"
    hd_ = deg2rad(hd_)-pi; % range(-pi,+pi)
elseif deg_or_rad == "deg"
    hd_ = hd_;
end

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

% break down 'refLoc2'
rlX2 = refLoc2(1,1);
rlY2 = refLoc2(1,2);

% compute egocentric angle to reference loc
% which one is correct?
% egoAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);

if deg_or_rad == "deg"
    egoAng = rem(atan2d(midY-rlY, midX-rlX)+180, 360);
elseif deg_or_rad == "rad"
    egoAng = rem(atan2d(midY-rlY, midX-rlX)+180, 360);
    egoAng = deg2rad(egoAng)-pi; % range(-pi,+pi)
end

% same, but for ref2
egoAng2 = rem(atan2d(midY-rlY2, midX-rlX2)+180, 360);
egoAng2 = deg2rad(egoAng2)-pi;


% find time indices when cell spikes
idx = knnsearch(t, SpikeTimes_);
spkX = x1(idx); spkY = y2(idx);
spk_egoAng = egoAng(idx);
spk_egoAng2 = egoAng2(idx); % for ref #2


%% compute tuning curves
nBins = 20;

if deg_or_rad == "rad"
    angEdges = linspace(-pi,pi,nBins);
elseif deg_or_rad == "deg"
    angEdges = linspace(0,360,nBins);
end
    

%% egoAng
[spkEgoAngMap, mapAxis_EgoAng] = histcounts(spk_egoAng,angEdges);
[allEgoAngMap] = histcounts(egoAng,angEdges);

for i = 1:length(mapAxis_EgoAng)
    if i+1 <= length(mapAxis_EgoAng)
        binCtrs_egoAng(i) = ((mapAxis_EgoAng(i+1)-mapAxis_EgoAng(i))/2)+mapAxis_EgoAng(i);
    end
end
tcVals_egoAng = spkEgoAngMap./(allEgoAngMap*sampleRate + eps); 

% smooth tuning curve values
tcVals_egoAng = imgaussfilt(tcVals_egoAng, 2, 'Padding', 'circular');


%% for reference point #2
[spkEgoAngMap2, mapAxis_EgoAng2] = histcounts(spk_egoAng2,angEdges);
[allEgoAngMap2] = histcounts(egoAng2,angEdges);

for i = 1:length(mapAxis_EgoAng2)
    if i+1 <= length(mapAxis_EgoAng2)
        binCtrs_egoAng2(i) = ((mapAxis_EgoAng2(i+1)-mapAxis_EgoAng2(i))/2)+mapAxis_EgoAng2(i);
    end
end
tcVals_egoAng2 = spkEgoAngMap2./(allEgoAngMap2*sampleRate + eps); 

% smooth tuning curve values
tcVals_egoAng2 = imgaussfilt(tcVals_egoAng2, 2, 'Padding', 'circular');


%% plot
if doPlot == "True"
    plot(binCtrs_egoAng, tcVals_egoAng, 'Color', 'k', 'LineWidth', 1.10)
    hold on
    plot(binCtrs_egoAng, tcVals_egoAng2, 'LineStyle', ':', 'Color', 'r', 'LineWidth', 1.10)
    % legend('center','hw', 'Location','northeastoutside', 'orientation', 'vertical')
    title("Egocentric Angle")
    ylabel("fr (Hz)")
    box off
    
    % decide how to label plot
    if deg_or_rad == "rad"
        xlim([-pi pi])
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
        xlabel("angle (rad)")
    elseif deg_or_rad == "deg"
        xlim([0 360])
        xticks([0 90 180 270 360])
        xlabel("angle (deg)")
    end
    
end
end

