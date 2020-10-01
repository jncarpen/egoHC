function [refVec] = generate_reference_pnts(pos_, extend_arena, nBins)
%GENERATE_REFERENCE_PNTS 
%   Generate a grid of reference points.

%   INPUTS
%   pos_cm:             position matrix [t x1 y1 x2 y2], with units in centimeters.
%
%   extend_arena:       'True' means you want to also get values from OUTSIDE of
%                       the arena
%   OUTPUTS 
%   refVec:             vector of reference points in the form [x y]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse position vector
    t = pos_(:,1);
    x = pos_(:,2);
    y = pos_(:,3);
    
    % Set default bin size to 10
%     nBins = 10;
    sampleRate = mode(diff(t));
    
    %% Bin spatial positions
    
    if extend_arena == "False"
         % Compute spatial occupancy and indices for X and Y bins (binX & binY)
        [~,xEdges,yEdges,~,~] = histcounts2(x,y,nBins);
    
    elseif extend_arena == "True"
        % Compute spatial occupancy and indices for extended X and Y bins
        scaleFac = .2; % how much to scale up arena (proportion)
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
end

