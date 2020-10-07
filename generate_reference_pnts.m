function [refVec, outside_circ] = generate_reference_pnts(position, extend_arena, nBins)
%GENERATE_REFERENCE_PNTS 
%   Generate a grid of reference points.

%   INPUTS
%   pos_cm:             position matrix [t x1 y1 x2 y2], with units in centimeters.
%
%   extend_arena:       'True' means you want to also get values from OUTSIDE of
%                       the arena
%   OUTPUTS 
%   refVec:             vector of reference points in the form [x y]
%   outside_circ:       points that are OUTSIDE of the radius you're
%                       interested in.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse position vector
    t = position(:,1);
    x = position(:,2);
    y = position(:,3);
    sampleRate = mode(diff(t));
    
    %% Bin spatial positions
    
    if extend_arena == "False"
         % Compute spatial occupancy and indices for X and Y bins (binX & binY)
        [~,xEdges,yEdges,~,~] = histcounts2(x,y,nBins);
    
    elseif extend_arena == "True"
        % Compute spatial occupancy and indices for extended X and Y bins
        scaleFac = 1.25; % how much to scale up arena (proportion)
        xSize = nanmax(x)-nanmin(x); ySize = nanmax(y)-nanmin(y);
        extended_x_edges = linspace(nanmin(x)-xSize*scaleFac, nanmax(x)+xSize*scaleFac, nBins+1);
        extended_y_edges = linspace(nanmin(y)-ySize*scaleFac, nanmax(y)+ySize*scaleFac, nBins+1);
        [~,xEdges,yEdges,~,~] = histcounts2(x,y,extended_x_edges,extended_y_edges);
    end

    
    %% Compute bin centers

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

    % make array of bin centers
    count = 1;
    for xx = 1:length(xCenter)
        for yy = 1:length(yCenter)
            refVec(count,1:2) = [xCenter(xx), yCenter(yy)];
            count = count+1;
        end
    end
    
    
    if extend_arena == "True"
        
        % define the circle of points we're interested in
        centerX = (nanmax(x)-nanmin(x))/2;
        centerY = (nanmax(y)-nanmin(y))/2;
        centers = [centerX, centerY];
        radius = 300;

        % compute the distance from each point from (centerX, centerY)
        refVec_circ = [];
        rvX = refVec(:,1); rvY = refVec(:,2);
        d = sqrt((centerX-rvX).^2 + (centerY-rvY).^2);
        in_idx = find(radius>=d); out_idx = find(radius<=d);
        outside_circ(:,1)= rvX(out_idx); outside_circ(:,2)= rvY(out_idx);
        refVec_circ(:,1) = rvX(in_idx); refVec_circ(:,2) = rvY(in_idx);
        refVec = refVec_circ;
        
    elseif extend_arena == "False"
        
        % shitty line of code cuz not parsing inputs well
        outside_circ = refVec;
    end
end

%% plot the reference points
% this can be used to check the reference points
% figure
% skip=5;
% for loc = 1:skip:length(refVec_circ)
%     plot(refVec_circ(loc,1), refVec_circ(loc,2), '.')
%     hold on
% end
% viscircles(centers,radius)
% plot(x,y)
% hold off


%% other scratch code:
        
%         rvX(radius>=d)=NaN; rvY(radius>=d)=NaN;
%         refVec = [rvX, rvY];


