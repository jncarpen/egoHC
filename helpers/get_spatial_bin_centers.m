function [binCtrs] = get_spatial_bin_centers(position, nBins)
%GET_SPATIAL_BIN_CENTERS Summary of this function goes here
%   INPUTS
%   'position'      [t x y x2 y2]
%   'nBins'         sqrt(# bins) to divide the spatial arena into. So if
%                   you want 100 bins, the input you want to pass is '10'.
%                   These obviously should be scalar values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse position input
t = position(:,1);
x = position(:,2);
y = position(:,3);

% Compute spatial occupancy and indices for X and Y bins (binX & binY)
[~,xEdges,yEdges,~,~] = histcounts2(x,y,nBins);

% compute spatial bin *centers* (x and y for quiverPlot)
for i = 1:length(xEdges)
    if i+1 <= length(xEdges)
        xCenter(i) = ((xEdges(i+1)-xEdges(i))/2)+xEdges(i);
    end
end

for i = 1:length(yEdges)
    if i+1 <= length(yEdges)
        yCenter(i) = ((yEdges(i+1)-yEdges(i))/2)+yEdges(i);
    end
end

count = 1;
for xx = 1:length(xCenter)
    for yy = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(xx), yCenter(yy)];
        count = count+1;
    end
end

% reshape the binCenters
binCtrs=[];
number = 1;
for ii = 1:10
    binCtrs= [binCtrs; binCenters(ii:10:end,:)];
    number = number + 1;
end

end



% more about reshaping the bin centers:
% the first point is the top left corner --> across to the
% top right corner and the second row starts at the left corner
% kind of like if you were reading

