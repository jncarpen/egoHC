function [HD_TC, ALLO_TC, EGO_TC, HD_ST, ALLO_ST, EGO_ST] = TC_stats_2DBins(position, head_direction, SpkTrn, refLoc, fileCount)
%TC_SCORING Summary of this function goes here
%   Get tuning curve statistics for each 2D bin for:
%   (1) HD tuning curve
%   (2) allocentric bearing tuning curve
%   (3) egocentric bearing tuning curve
%   There are 10x10 total spatial bins and a tuning curve is computed for
%   each.
%   Jo Carpenter, October 5, 2020.

%% Set things up

% parse position vector
t = position(:,1);
x1 = position(:,2); y1 = position(:,3);
x2 = position(:,4); y2 = position(:,5);

% compute 2D spatial bins
nBins = 10; 
sampleRate = mode(diff(t));

% Compute spatial occupancy and indices for X and Y bins (binX & binY)
% Don't remember where I used this next??
[spatialOcc,xEdges,yEdges,binX,binY] = histcounts2(x1,y1,nBins);

%% calculate egocentric bearing (for all timepoints)

% find midpoint between two LEDs
midX = (x1+x2)/2; midY = (y1+y2)/2;

% break down refLoc
rlX = refLoc(1,1); rlY = refLoc(1,2);

% find allocentric + egocentric 'bearing' at each timepoint (AB/EB)
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
egoAng = alloAng - head_direction;
% correct for negative angles (egoAng)
neg_idx = find(egoAng<0);
egoAng(neg_idx) = egoAng(neg_idx)+360;

%% compute spatial bin centers 
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

% make a vector of the bin centers
count = 1;
for xx = 1:length(xCenter)
    for yy = 1:length(yCenter)
        ctrLocs(count,1:2) = [xCenter(xx), yCenter(yy)];
        count = count+1;
    end
end

%% Compute a tuning curves for each 2D spatial bin
numBins_HD = 40; % 9 degree bins
binWidth_deg = 360/numBins_HD;
angBins = linspace(0,360,numBins_HD);

% set up all the cell arrays youre gonna need
HD_TC = cell(10,10); EGO_TC = cell(10,10); ALLO_TC = cell(10,10);
HD_ST = cell(10,10); EGO_ST = cell(10,10); ALLO_ST = cell(10,10);

% iterate through every 2D spatial bin 
count = 1;
for xx = 1:nBins
    for yy = 1:nBins
        indices = find(xx == binX & yy == binY);
        timeInBin = length(indices)*sampleRate; % occupancy (s)

        % calculate values for current 2D spatial bin
        spikes_here = SpkTrn(indices); hd_here = head_direction(indices); 
        allo_here = alloAng(indices); ego_here = egoAng(indices);

        % find HD, AB, EB at time of spikes 
        spikeInds = find(spikes_here > 0); % since its smoothed
        angSpk = hd_here(spikeInds); 
        alloSpk = allo_here(spikeInds); egoSpk = ego_here(spikeInds);

        % compute normal firing rate map [r(x,y)] for this 2D spatial bin
        rateMap_HD(xx,yy) = sum(SpkTrn(indices))./(timeInBin);
        
        % compute tuning curves
        hd_tc = analyses.turningCurve(angSpk, hd_here, sampleRate, 'smooth', 1, 'binWidth', binWidth_deg);
        allo_tc = analyses.turningCurve(alloSpk, allo_here, sampleRate, 'smooth', 1, 'binWidth', binWidth_deg);
        ego_tc = analyses.turningCurve(egoSpk, ego_here, sampleRate, 'smooth', 1, 'binWidth', binWidth_deg);
        
        % compute tuning curve statistics
        hd_tcStat = analyses.tcStatistics(hd_tc, binWidth_deg, 95);
        allo_tcStat = analyses.tcStatistics(allo_tc, binWidth_deg, 95);
        ego_tcStat = analyses.tcStatistics(ego_tc, binWidth_deg, 95);
        
        
        % put everything into cell arrays (10x10)
        HD_TC{xx,yy} = hd_tc; ALLO_TC{xx,yy} = allo_tc; EGO_TC{xx,yy} = ego_tc;
        HD_ST{xx,yy} = hd_tcStat; ALLO_ST{xx,yy} = allo_tcStat; EGO_ST{xx,yy} = ego_tcStat;
        
        % get values for peak direction
        hd_PD(count)=hd_tcStat.peakDirection;
        allo_PD(count)=allo_tcStat.peakDirection;
        ego_PD(count)=ego_tcStat.peakDirection;
        
        % get values for mean direction
        hd_MD(count)=hd_tcStat.mean;
        allo_MD(count)=allo_tcStat.mean;
        ego_MD(count)=ego_tcStat.mean;
        
        % get values for MVL
        hd_MVL(count)=hd_tcStat.r;
        allo_MVL(count)=allo_tcStat.r;
        ego_MVL(count)=ego_tcStat.r;
        count = count+1;
    end
end

%% plot the results
theta = ego_PD; theta2 = hd_PD;
quiver_handle = plot_quiver_stacked(theta, theta2, ctrLocs); % quiver function is already holding
h1 = plot(rlX, rlY, 'o', 'MarkerSize', 12);
set(h1, 'markerfacecolor', 'red');

% save the figure
% filename = strcat('D:\egoAnalysis\ego_hd_peakRate\', 'refX', sprintf('%.f', rlX), '_refY', sprintf('%.f', rlY), '.png');
filename = strcat('D:\egoAnalysis\ego_hd_peakRate2\', sprintf('%.f', fileCount), '.png');
saveas(quiver_handle, filename);
end

