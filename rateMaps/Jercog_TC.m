function Jercog_TC(pos_, hd_, SpikeTrain)
%HDMOD: Jercog et al. analysis knockoff

    %   INPUT:
    %   pos:                Position samples, [t x y]
    %   hd:                 Vector of head direction values
    %   SpikeTrain:         Spike train for cell

    %   OUTPUT:
    %   A bunch of stuff right now, will be refined once function is
    %   complete.
    
    %   Jordan Carpenter 2020

%% Set everything up
    
    % Parse position vector
    t = pos_(:,1);
    x = pos_(:,2);
    y = pos_(:,3);
    
    % convert HD from degrees to radians (-pi:+pi)
    hd_rad = deg2rad(hd_)-pi;
    
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
            scatter(xCenter(xx), yCenter(yy)) % plot vector locations (to
%             check)
            count = count+1;
            hold on
        end
    end
    
       

%%  Compute tuning curves for each 2D spatial bin

    % compute bins (-pi:+pi)
    numBins_HD = 10;
    angBins = linspace(-pi,pi,numBins_HD);
    
    % initialize cell array (tuning curve values)
    hd_tc = cell(10,10);
    
    for xx = 1:nBins
        for yy = 1:nBins
            indices = find(xx == binX & yy == binY);
            timeInBin = length(indices)*sampleRate; % amount of time spent in this 2D spatial bin (occupancy)
            
            % calculate values for current 2D spatial bin
            hd_here = hd_rad(indices); % head direction values in this spatial bin
%             time_here = t(indices); % time values in this spatial bin
            spikes_here = SpikeTrain(indices); % spiketrain in this spatial bin
            
            % find head-direction at times of spikes (angSpk)
            
            spikeInds = find(spikes_here >= 1);
            angSpk = hd_here(spikeInds);
            
%             count = 1; 
%             for index = 1:length(spikes_here)
%                 if spikes_here(index) >= 1
%                     angSpk(count)= hd_here(index);
%                 end
%             end
            
            % compute firing rate map
            rateMap(xx,yy) = sum(SpikeTrain(indices))/(timeInBin); % normal rate map [r(x,y)]

            % compute occupancy for each HD bin
            histAng = histcounts(hd_rad(indices), angBins);
            spkPerAng = histcounts(angSpk,angBins);
            hd_tc{xx,yy} = spkPerAng./histAng * sampleRate;
        end
    end
    
    % compute angular bin centers
    for i = 1:length(angBins)
        if i+1 <= length(angBins)
            angBins_ctr(i) = ((angBins(i+1)-angBins(i))/2)+angBins(i);
        end
    end
           
    
    %% plot stuff!
    
    peakRate = max(max(rateMap)); % define peak FR (across ALL 2D spatial bins)
    
    % compute maximum tuning curve value (in the cartesian plane)
    for xx = 1:10
        for yy = 1:10
            tcHere = hd_tc{xx,yy};
            tcHere(isnan(tcHere)) = 0; % set NaNs to zero
            [angBins_cart,hd_tc_cart] = pol2cart(angBins_ctr, tcHere);
            maxTC(xx,yy) = max(max(hd_tc_cart));
            maxBins(xx,yy) = max(max(angBins_cart));
            
            minTC(xx,yy) = min(min(hd_tc_cart));
            minBins(xx,yy) = min(min(angBins_cart)); 
        end
    end
    maxTC_cart = max(max(maxTC)); 
    maxBin_cart = max(max(maxBins));
    
    minTC_cart = min(min(minTC)); 
    minBin_cart = min(max(minBins));
    
    
    
    
    % f = figure;
%     set(f,'Position',[40 40 600 600])
    count = 1;
    for xx = 5:10
        for yy = 5:10
            % convert polar to cartesian coordinates
            tcHere = hd_tc{xx,yy};
            tcHere(isnan(tcHere)) = 0; % set NaNs to zero
            [angBins_cart,hd_tc_cart] = pol2cart(angBins_ctr, tcHere);
            
            % compute HD tuning curve for spatial bins
            sub = subplot(10,10,count);
            imagesc(rateMap(xx,yy))
            caxis([0, peakRate]) % normalize color axis
            set(gca,'YDir','normal')
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            p = get(gca, 'Position');
            h = axes('Parent', gcf, 'Position', [p(1) p(2) p(3) p(4)]);
            set(h,'Color','none')
            hold on
            plot(angBins_cart, hd_tc_cart, 'k', 'LineWidth', 1.25)
            ylim([minTC_cart, maxTC_cart])
            xlim([minBin_cart, maxBin_cart])
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            count = count + 1;
            
        end
    end
end
