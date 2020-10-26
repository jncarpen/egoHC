function [tcVals_egoAng, binCtrs_egoAng] = egoBearing(position, ST, refLoc, refLoc2, doPlot, deg_or_rad)
%EGOBEARING Compute tuning curve for egocentric bearing

% calculate head direction
[head_direction] = get_hd(position);

% decide whether to convert to radians or keep in degrees
if deg_or_rad == "rad"
    head_direction = deg2rad(head_direction)-pi; % range(-pi,+pi)
elseif deg_or_rad == "deg"
    head_direction = head_direction;
end

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

% break down 'refLoc2'
rlX2 = refLoc2(1,1);
rlY2 = refLoc2(1,2);

if deg_or_rad == "deg"
    alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
    egoAng = alloAng - head_direction;
    % correct for negative angles if any are found
    neg_idx = find(egoAng<0);
    egoAng(neg_idx) = egoAng(neg_idx)+360;
    
%     egoAng = mod(egoAng, 360);
%     egoAng = rem(egoAng, 360);
    
elseif deg_or_rad == "rad"
    alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
    egoAng = deg2rad(alloAng-head_direction)-pi; % or is it 2pi?
end

% same, but for ref2
if deg_or_rad == "deg"
    % this is the way i computed it in the python script
    alloAng2 = rem(atan2d(rlY2-midY, rlX2-midX)+180, 360);
    egoAng2 = alloAng2 - head_direction;
    % correct for negative angles if any are found
    neg_idx2 = find(egoAng2<0);
    egoAng2(neg_idx2) = egoAng2(neg_idx2)+360;
    
elseif deg_or_rad == "rad"
    alloAng2 = rem(atan2d(rlY2-midY, rlX2-midX)+180, 360);
    egoAng2 = deg2rad(alloAng2-head_direction)-pi;
end

%% compute tuning curves
nBins = 40; % 9 degree bins

if deg_or_rad == "rad"
    angEdges = linspace(-pi,pi,nBins);
elseif deg_or_rad == "deg"
    angEdges = linspace(0,360,nBins);
end

% find time indices when cell spikes
idx = knnsearch(t, ST);
spkX = x1(idx); spkY = y2(idx);
spk_egoAng = egoAng(idx);
spk_egoAng2 = egoAng2(idx);

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
    figure
    plot(binCtrs_egoAng, tcVals_egoAng, 'Color', 'k', 'LineWidth', 1.10)
    hold on
    plot(binCtrs_egoAng2, tcVals_egoAng2, 'LineStyle', ':', 'Color', 'r', 'LineWidth', 1.10)
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

%% SCRATCH CODE
% grab stuff
% if sessNum == "False"
%     pos_ = pos_cm;
%     SpikeTimes_ = SpikeTimes;
%     hd_ = hd;
% else
% end

% compute egocentric angle to reference loc
% which one is correct?
% egoAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);

% this is the way I was originally computing ego angle.
% egoAng = rem(atan2d(midY-rlY, midX-rlX)+180, 360);
