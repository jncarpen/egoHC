function occ = occupancyMap(P, nBins)
%OCCUPANCYMAP Summary of this function goes here
%   Detailed explanation goes here
% plot occ and map
posx = P(:,2);
posy = P(:,3);
tpf = mode(diff(P(:,1))); % time per frame

% get min/max values
xMin = nanmin(P(:,2));
xMax = nanmax(P(:,2));
yMin = nanmin(P(:,3));
yMax = nanmax(P(:,3));

maxEdges = max(xMax, yMax);
minEdges = min(xMin, yMin);

% bin position data
% nBins = 10;
[occ, xEdges, yEdges, binX, binY] = histcounts2(posx,posy,nBins);

% for rr = 1:nBins
%     for cc = 1:nBins
%         % find frames in which animal occupied this spatial bin
%         occ(rr,cc) = length(find(rr == binY & cc == binX));
%     end
% end
% 
% figure
% imagesc(occ)
% pbaspect([1 1 1])
% title("Spatial Occ")
% xlabel("X")
% ylabel("Y")
% colormap('jet')
% colorbar


return

end

