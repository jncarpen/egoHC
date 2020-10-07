function [countMap, occMap, rateMap] = FRmap(pos, SpikeTrain)
%FRMAP Generate and plot firing rate map

%   DESCRIPTION:
%   2D location is divided into 2x2cm bins. The number of spikes fired when
%   the animal occupies each bin is calculated. This is divided by the amount
%   of time (in seconds) that the animal occupies each bin. The resulting
%   heatmap is smoothed with a Gaussian kernel with standard deviation = 2
%   bins.

%   INPUT:
%   pos:            [t x y x2 y2]
%   SpikeTimes:     Binned spikes; length(pos) == length(SpikeTimes)

%   OUTPUT:

%   Jordan Carpenter, 2020.

%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%
    % Make sure BNT is on the path (for smooth function)
    addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001\bnt"));
    
    % Parse position vector
    t = pos(:,1);
    x = pos(:,2);
    y = pos(:,3);
    
    % Set default bin size to 10
    nBins = 50;
    sampleRate = mode(diff(t));
    
    % Compute spatial occupancy and indices for X and Y bins (binX & binY)
    [spatialOcc,xEdges,yEdges,binX,binY] = histcounts2(x,y,nBins);
    
     for xx = 1:nBins
        for yy = 1:nBins
            indices = find(xx == binX & yy == binY);
            countMap(xx,yy) = sum(SpikeTrain(indices));
            occMap(xx,yy) = (length(indices)*sampleRate);
%             rateMap(xx,yy) = sum(SpikeTrain(indices))/(length(indices)*sampleRate);
        end
     end
     
     % Smooth occupacy map a bit to deal with outliers...
     % USAGE: smoothed = general.smooth(data, [verticalSD, horizSD]);
     smoothOCC = general.smooth(occMap, [2,2]);
     
     % divide CM by occu to get rate map
     rateMap = countMap./smoothOCC;
     
     
     

