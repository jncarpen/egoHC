function [R, weightSum,binX, vecLocs, rateMap, rateMap_cond] = hdMod(pos_, hd_, SpikeTrain_, sessNum, unitNum)
%HDMOD: Jercog et al. analysis knockoff :p

    %   INPUT:
    %   pos:                Position samples, [t x y] for a single animal
    %   hd:                 Vector of head direction values (in *degrees*) for a single animal
    %   SpikeTrain:         Spike train for cell for a single animal
    %   sessNum:            Session number you want to process.
    %   unitNum:            Unit number you want to process.

    %   OUTPUT:
    %   A bunch of stuff right now, will be refined once function is
    %   complete.
    
    %   Modified by J. Carpenter 2020

%% Set everything up

    % Grab correct session/unit
    pos_ = pos{1,sessNum};
    hd_ = deg2rad(hd{1,sessNum})-pi; % convert to radians (-pi:+pi)
    SpikeTrain_ = SpikeTrain{1,sessNum}{1,unitNum};
    
    % Parse position vector
    t = pos_(:,1);
    x = pos_(:,2);
    y = pos_(:,3);
    
    % Set default bin size to 10
    nBins = 10;
    sampleRate = mode(diff(t));
    
    % Compute spatial occupancy and indices for X and Y bins (binX & binY)
    [spatialOcc,xEdges,yEdges,binX,binY] = histcounts2(x,y,nBins);
    
    % 'binX' and 'binY': Arrays the same size as 'x' and 'y' whose elements are the
    % bin indices for the corresponding elements in X. I am going to use
    % these to compute the spatial occupancy for each index by counting the
    % number of '3's there are, for example, and dividing that by the
    % sampleRate that was calculated. This will give me the total amount of
    % time the animal spent in a particular 2D spatial bin.
    
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
    
    % make scatter plot of grids
    % this was just for me to check
    count = 1;
    for xx = 1:length(xCenter)
        for yy = 1:length(yCenter)
            vecLocs(count,1:2) = [xCenter(xx), yCenter(yy)];
%             scatter(xCenter(xx), yCenter(yy)) % plot vector locations (to
%             check)
            count = count+1;
            hold on
        end
    end
    
       

%%  determine preferred HD vector
    angBins = 10; % define number of angular HD bins
    binEdgesHD = linspace(-pi,pi,angBins+1); % calculate bin edges (-pi:+pi)
    
    
    % compute angular bin centers (these are the 10 unit vectors for each
    % spatial bin pointing in the directions defined by the centers of the
    % 10 spatial bins).
    for i = 1:length(binEdgesHD)
        if i+1 <= length(binEdgesHD)
            HDBins_ctr(i) = ((binEdgesHD(i+1)-binEdgesHD(i))/2)+binEdgesHD(i);
        end
    end
           
    
    for xx = 1:nBins
        for yy = 1:nBins
            indices = find(xx == binX & yy == binY);
            
            % compute mean HD in each spatial bin (vector *direction*)
            HDmean(xx, yy) = circ_mean(hd_(indices));
            
            % meanFR(xx, yy) = mean(cellSpk(indices)); % don't remember why
            % this is here...
            
            % compute maps
            timeInBin = length(indices)*sampleRate; % amount of time spent in this 2D spatial bin (occupancy)
            rateMap(xx,yy) = sum(SpikeTrain_(indices))/(timeInBin); % normal rate map [r(x,y)]
%             occuMap(xx,yy) = length(indices);
%             countMap(xx,yy) = sum(SpikeTrain(indices));
            
            count = 1; % initialize count
            for H = binEdgesHD % loop through HD bin edges (that were defined in the HDbins vector)
                indicesHD = find(xx == binX & yy == binY & hd_>H & hd_<=H+(2*pi/10));
                timeInBin_HD = length(indicesHD)*sampleRate; % occupancy(spatial+HD)
                rateMap_cond(xx, yy, count) = sum(SpikeTrain_(indicesHD))/timeInBin_HD;
                count = count+1; % add 1 to count
            end
        end
    end
    
    % compute directional modulation ratio [R] (independent of spatial
    % modulation): R(x,y,H) = r(x,y,H)/r(x,y)
    % NOTE: *ADD* restriction to spatial bins r(x,y)>0.5Hz to avoid dividing by 0.
    R = rateMap_cond./rateMap;
    R(isnan(R))=0; % set NaN values to zero (this might be a mistake)
    
    % plot this (just to visualize for now)
    for bin = 1:10
        subplot(5,2,bin)
        imagesc(R(:,:,bin))
        set(gca,'YDir','normal') 
        colormap(jet)
        colorbar
        caxis([1, 30])
    end
    
    % compute weighted sum
    weightSum = R.* binEdgesHD;
    
    % sum over H (note: weightSum(x,y,H)) to determine preferred
    % HD vector for each spatial bin
    
    for xx = 1:nBins
        for yy = 1:nBins
            prefHDvect(xx,yy) = nansum(squeeze(weightSum(xx,yy,:)))/nBins;
        end
    end
    
    %% plot stuff!
    % this part will be removed later and replaced by just the plot we want
    % to make \\
    
    
    % plot rate map
    % Need to find a good filter \\
    % imagesc(smoothdata(rateMap, 'gaussian',10));
    % imagesc(conv2(countMap, ones(11,1)/11, 'same')); % moving avg
    
    % sample quiver components (change later)
    u = cos(vecLocs(:,1)).*vecLocs(:,2);
    v = cos(vecLocs(:,1)).*vecLocs(:,2);
    
    figure
%     imagesc(rateMap); % plot ratemap
%     hold on
    quiver(vecLocs(:,1), vecLocs(:,2), u, v)
    
    set(gca,'YDir','normal') % inverts the y-axis direction (to normal)
    title("Firing Rate Map")
    xlabel("XPos")
    ylabel("YPos")
    

% EXTRA PLOTS (maybe move to another function)    
%     subplot(2,2,2)
%     imagesc(countMap);
%     set(gca,'YDir','normal')
%     title("Spike Counts")
%     xlabel("XPos")
%     ylabel("YPos")
%     
%     subplot(2,2,3)
%     imagesc(occuMap);
%     set(gca,'YDir','normal')
%     title("Spatial Occupancy")
%     xlabel("XPos")
%     ylabel("YPos")


%% plot unit vectors
% Plot a circle.
angles = linspace(-pi, pi, 10); % 720 is the total number of points
radius = 1;
xCenter = 0;
yCenter = 0;
x = radius * cos(angles) + xCenter; 
y = radius * sin(angles) + yCenter;
% plot center
plot(xCenter, yCenter, 'k+', 'LineWidth', 2, 'MarkerSize', 5);
hold on
grid on;
axis equal;
xlabel('X', 'FontSize', 10);
ylabel('Y', 'FontSize', 10);
for lineNum = 1:10
    A = [xCenter, x(lineNum)];
    B = [yCenter, y(lineNum)];
    plot(A, B, 'k','LineWidth', 1.25);
end
    
end