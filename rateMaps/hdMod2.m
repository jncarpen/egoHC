function hdMod2(pos, hd, SpikeTrain, sessNum, unitNum)
%HDMOD: Jercog et al. analysis knockoff :p (V2)

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
    %   hdMod V2

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
    count = 1;
    for xx = 1:length(xCenter)
        for yy = 1:length(yCenter)
            vecLocs(count,1:2) = [xCenter(xx), yCenter(yy)];
            scatter(xCenter(xx), yCenter(yy)) % plot vector locations (to
%             check)
            hold on
            count = count+1;
        end
    end
    
    %% Determine preferred HD vector
    angBins = 10; % define number of angular HD bins
    binEdgesHD = linspace(-pi,pi,angBins+1); % calculate bin edges (-pi:+pi)
    
   % Compute binEdgesHD centers
   for i = 1:length(binEdgesHD)
        if i+1 <= length(binEdgesHD)
            HDBins_ctr(i) = ((binEdgesHD(i+1)-binEdgesHD(i))/2)+binEdgesHD(i);
        end
    end
    
    
    % loop through each 2D spatial bin
    for xx = 1:nBins
        for yy = 1:nBins
            
            % find indices in which animal occupied this 2D spatial bin
            indices = find(xx == binX & yy == binY);
            
            % compute mean HD in this 2D spatial bin
            HDmean(xx, yy) = circ_mean(hd_(indices));
            
            % amount of time spent in 2D spatial bin
            timeInBin = length(indices)*sampleRate;
            
            % compute spatial rate map- r(x,y)
            rateMap(xx,yy) = sum(SpikeTrain_(indices))/timeInBin;
            
            % iterate through HD bins
            count = 1;
            for H = 1:length(binEdgesHD)-1 % loop through HD bin edges (that were defined in the HDbins vector)
                
                % find indices in which animal occupied this HD bin
                indicesHD = find(xx == binX & yy == binY & hd_>binEdgesHD(H) & hd_<binEdgesHD(H+1));
                
                % amount of time spent in HD bin
                timeInBin_HD = length(indicesHD)*sampleRate; % occupancy(spatial+HD)
                
                % compute conditional rate map- r(x,y,H)
                rateMap_cond(xx, yy, count) = sum(SpikeTrain_(indicesHD))/timeInBin_HD;
                
                % directional modulation independent of spatial modulation
                R(xx, yy, count) = rateMap_cond(xx,yy,count)./rateMap(xx,yy);
                
                % weight unit vector (HDBins_ctr(H)) by R(x,y,H)
                weightVec(xx, yy, count) = HDBins_ctr(count).*R(xx, yy, count);
                
                count = count+1;
            end 
        end
    end
    
 % set NaNs to zero
 R(isnan(R))=0;
 weightVec(isnan(weightVec))=0; 
 
 % calculate preferred HD vector
  for xx = 1:nBins
      for yy = 1:nBins
          prefVec(xx,yy) = sum(weightVec(xx,yy,:)); 
          prefR(xx,yy) = sum(R(xx,yy,:));
      end
  end  
  
  
    
end
    
    
    
    