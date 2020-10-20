function [tcVals_egoAng] = egoBearing_heading(pos_, STimes, refLoc, refLoc2, doPlot, deg_or_rad, method)
%EGOBEARING Compute egocentric bearing tuning curves (but with respect to
%the *heading* direction using one of three measures.

% compute moving direction
[md_1, md_2, md_3] = get_moving_direction(pos_);

% choose which method to use
switch method
    case "1"
        heading_dir = md_1';
    case "2"
        heading_dir = md_2';
    case "3"
        heading_dir = md_3';
end
        

% decide whether to convert to radians or keep in degrees
if deg_or_rad == "rad"
    heading_dir = deg2rad(heading_dir)-pi; % range(-pi,+pi)
elseif deg_or_rad == "deg"
    heading_dir = heading_dir;
end

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

if deg_or_rad == "deg"
    % this is the way i computed it in the python script
    alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
    egoAng = alloAng - heading_dir;
    % correct for negative angles if any are found
    neg_idx = find(egoAng<0);
    egoAng(neg_idx) = egoAng(neg_idx)+360;
    
elseif deg_or_rad == "rad"
    alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
    egoAng = deg2rad(alloAng-heading_dir)-pi;
end

% same, but for ref2
if deg_or_rad == "deg"
    % this is the way i computed it in the python script
    alloAng2 = rem(atan2d(rlY2-midY, rlX2-midX)+180, 360);
    egoAng2 = alloAng2 - heading_dir;
    % correct for negative angles if any are found
    neg_idx2 = find(egoAng2<0);
    egoAng2(neg_idx2) = egoAng2(neg_idx2)+360;
    
elseif deg_or_rad == "rad"
    alloAng2 = rem(atan2d(rlY2-midY, rlX2-midX)+180, 360);
    egoAng2 = deg2rad(alloAng2-heading_dir)-pi;
end


% find time indices when cell spikes
idx = knnsearch(t, STimes);
% spkX = x1(idx); spkY = y2(idx);
spk_egoAng = egoAng(idx);
spk_egoAng2 = egoAng2(idx); % for ref #2


%% compute tuning curves
nBins = 40; % 9 degree bins

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
%     figure
    plot(binCtrs_egoAng, tcVals_egoAng, 'Color', 'k', 'LineWidth', 1.10)
    hold on
    plot(binCtrs_egoAng2, tcVals_egoAng2, 'LineStyle', ':', 'Color', 'r', 'LineWidth', 1.10)
    % legend('center','hw', 'Location','northeastoutside', 'orientation', 'vertical')
    title("Egocentric Angle (Heading Dir)")
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

