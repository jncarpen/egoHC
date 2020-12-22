function [initial, grid] = choose_initial_conditions_grid(P, howMany)

% break up position vector
x = P(:,2); y = P(:,3);

% divide the arena into 100 2D spatial bins
[~, xEdges, yEdges, ~, ~] = histcounts2(x,y,howMany);
yEdges = fliplr(yEdges); % flip y-vector

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

% count = 1;
for col = 1:length(xCenter)
    for row = 1:length(yCenter)
%         binCenters(count,1:2) = [xCenter(row), yCenter(col)];
        xmat(col,row) = xCenter(row);
        ymat(col,row) = yCenter(col);
        
        initial{col,row}.g = .5;
        initial{col,row}.thetaP = 0;
        initial{col,row}.xref = xCenter(row);
        initial{col,row}.yref = yCenter(col);

%         count = count+1;
    end
end

% save matrix stuff
grid.x = xmat;
grid.y = ymat;



end

