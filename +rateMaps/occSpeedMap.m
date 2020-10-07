function [occ, map] = occSpeedMap(pos, speed)
%OCCUPANCYMAP: create occupancy and speed maps

% plot occ and map
posx = pos(:,2);
posy = pos(:,3);

% DISCRETIZE X and Y

% bin position data
posBins = 20;
posEdges = linspace(0, 80, posBins+1);

% FOR X:
% compute position occupancy & pos indices
[posXOccupancy,~,posXInds] = histcounts(posx, posEdges);

% initialize zero matrix
posX_dis = zeros(length(posXInds), posBins);

for T=1:length(posXInds) % for every time point
    currentBin = posXInds(T); % pulls out value for HD bin at time_i
    posX_dis(T, currentBin+1) = 1; % sets value to 1
end


% FOR Y:

[posYOccupancy,~,posYInds] = histcounts(posy, posEdges);

% initialize zero matrix
currentBin = [];
posY_dis = zeros(length(posYInds), posBins);

for T=1:length(posYInds) % for every time point
    currentBin = posYInds(T); % pulls out value for posY
    posY_dis(T, currentBin+1) = 1; % sets value to 1
end

%% OCC
SF = 5; % smoothing factor
for i=1:100
    for j= 1:100
        f = find(abs(posX_dis-i)<SF & abs(posY_dis-j)<SF);
        occ(i,j) = length(f);
        map(i,j) = nanmedian(speed(f));
    end
end

figure 
imagesc(occ)
title("Occupancy")
xlabel("X")
ylabel("Y")
colorbar

% figure
% contourf(occ)
% title("Occupancy- Contour")
% xlabel("X-Coordinate")
% ylabel("Y-Coordinate")
% colorbar

figure
imagesc(map)
title("Speed Map")
xlabel("X")
ylabel("Y")
colorbar

% figure
% contourf(map)
% title("Speed Map- Contour")
% xlabel("X-Coordinate")
% ylabel("Y-Coordinate")
% colorbar

return
end

