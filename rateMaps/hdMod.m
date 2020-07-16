function [HDmean, rateMap, spkRate, weightSum, prefHDvect, occuMap,countMap] = hdMod(pos, head_angle, cellSpk)
%BINPOS: Jercog et al. analysis knockoff

    %   INPUT:
    %   pos             Position samples, [t x y]
    %   head_angle      Vector of head direction values
    %   cellSpk         Spike train for cell

    %   OUTPUT:
    %   A bunch of stuff right now, will be refined once function is
    %   complete.


%% Set everything up
    
    % Parse position vector
    t = pos(:,1);
    x = pos(:,2);
    y = pos(:,3);
    
    % Set default bin size to 10
    nBins = 10;
    sampleRate = mode(diff(t));
    
    % Compute spatial occupancy and indices for X and Y bins (binX & binY)
    [spatialOcc,xEdges,yEdges,binX,binY] = histcounts2(x,y,nBins);
    
    % compute spatial bin *centers* (for when we plot the HD-mod vectors)
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
       

%%  determine preferred HD vector

    HDbins = linspace(0,2*pi,10);
    for xx = 1:nBins
        for yy = 1:nBins
            indices = find(xx == binX & yy == binY);
            HDmean(xx, yy) = circ_mean(head_angle(indices)); % mean HD in each spatial bin
            % meanFR(xx, yy) = mean(cellSpk(indices));
            rateMap(xx,yy) = sum(cellSpk(indices))/(length(indices)*sampleRate/1000);
            occuMap(xx,yy) = length(indices);
            countMap(xx,yy) = sum(cellSpk(indices));
            count = 1;
            for H = HDbins
                indicesH = find(xx == binX & yy == binY & head_angle>H & head_angle<=H+(2*pi/10));
                spkRate(xx, yy, count) = sum(cellSpk(indicesH))/(length(indicesH)*sampleRate/1000);
                count = count+1;
            end
        end
    end
    
    % compute directional modulation ratio
    % NOTE: *ADD* restriction to spatial bins r(x,y)>0.5Hz to avoid dividing by 0.
    
    R = spkRate./rateMap;
    weightSum = R.* HDbins;
    
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
    
    figure
    
    % plot rate map
    subplot(2,2,1);
    % Need to find a good filter \\
    % imagesc(smoothdata(rateMap, 'gaussian',10));
    % imagesc(conv2(countMap, ones(11,1)/11, 'same')); % moving avg
    imagesc(rateMap); 
    title("Firing Rate Map")
    xlabel("XPos")
    ylabel("YPos")
    
    subplot(2,2,2)
    imagesc(countMap); 
    title("Spike Counts")
    xlabel("XPos")
    ylabel("YPos")
    
    subplot(2,2,3)
    imagesc(occuMap);
    title("Spatial Occupancy")
    xlabel("XPos")
    ylabel("YPos")
    
end