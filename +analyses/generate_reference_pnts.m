function [refVec, in_out_index] = generate_reference_pnts(position, extend_arena, nBins)
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

    
%% plot xedges
% plot(x,y)
% hold on
% for ll = 1:length(extended_x_edges)
%     xline(extended_x_edges(ll))
%     yline(extended_y_edges(ll))
% end
% hold off
    
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
    
    
    
    %% define points inside of the circle    
    
    % this vector will be filled with 1s and 0s; 
    % 1 means that the point is INSIDE the circle, 0 means its OUTSIDE
    in_out_index = zeros(length(refVec),1);
    
    if extend_arena == "True"
        
        % define the circle of points we're interested in
        centerX = (nanmax(x)-nanmin(x))/2;
        centerY = (nanmax(y)-nanmin(y))/2;
        centers = [centerX, centerY];
        radius = (abs(nanmin(x)-xSize*scaleFac) + abs(nanmax(x)+xSize*scaleFac))/2;
        

        % compute the distance from each point from (centerX, centerY)
%         refVec_circ = [];
        rvX = refVec(:,1); rvY = refVec(:,2);
        d = sqrt((centerX-rvX).^2 + (centerY-rvY).^2);
        in_idx = find(radius>=d); out_idx = find(radius<=d);
        in_out_index(in_idx) = 1; in_out_index(out_idx) = 0;
%         refVec_circ(:,1) = rvX(in_idx); refVec_circ(:,2) = rvY(in_idx);
%         refVec = refVec_circ;
        
    elseif extend_arena == "False"
        outside_circ = refVec; % shitty line of code cuz not parsing inputs well (sorry future self)
    end
end

% %% plot the reference points
% % this can be used to check the reference points
% figure
% skip=1;
% for loc = 1:skip:length(refVec)
%     plot(refVec(loc,1), refVec(loc,2), '.')
%     hold on
% end
% 
% viscircles(centers,radius)
% plot(x,y)
% hold off


%% other scratch code:
        
%         rvX(radius>=d)=NaN; rvY(radius>=d)=NaN;
%         refVec = [rvX, rvY];


