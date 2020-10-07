function [HD_TC, ALLO_TC, EGO_TC, HD_ST, ALLO_ST, EGO_ST] = TC_stats_2DBins(position, head_direction, SpkTrn, refLoc)
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
egoAng = abs(alloAng - head_direction); % is this correct?


%% Compute a tuning curves for each 2D spatial bin
numBins_HD = 40; % 9 degree bins
angBins = linspace(0,360,numBins_HD);

% set up all the cell arrays youre gonna need
HD_TC = cell(10,10); EGO_TC = cell(10,10); ALLO_TC = cell(10,10);
HD_ST = cell(10,10); EGO_ST = cell(10,10); ALLO_ST = cell(10,10);

% iterate through every 2D spatial bin 
for xx = 1:nBins
    for yy = 1:nBins
        indices = find(xx == binX & yy == binY);
        timeInBin = length(indices)*sampleRate; % occ (s)

        % calculate values for current 2D spatial bin
        spikes_here = SpkTrn(indices); hd_here = head_direction(indices); 
        allo_here = alloAng(indices); ego_here = egoAng(indices);

        % find HD, AB, EB at time of spikes 
        spikeInds = find(spikes_here >= 1);
        angSpk = hd_here(spikeInds); 
        alloSpk = allo_here(spikeInds); egoSpk = ego_here(spikeInds);

        % compute normal firing rate map [r(x,y)] for this 2D spatial bin
        rateMap_HD(xx,yy) = sum(SpkTrn(indices))./(timeInBin);
        
        % compute tuning curves
        hd_tc = analyses.turningCurve(angSpk, hd_here, sampleRate, 'smooth', 1, 'binWidth', 3);
        allo_tc = analyses.turningCurve(alloSpk, allo_here, sampleRate, 'smooth', 1, 'binWidth', 3);
        ego_tc = analyses.turningCurve(egoSpk, ego_here, sampleRate, 'smooth', 1, 'binWidth', 3);
        
        % compute tuning curve statistics
        hd_tcStat = analyses.tcStatistics(hd_tc, 3, 20);
        allo_tcStat = analyses.tcStatistics(allo_tc, 3, 20);
        ego_tcStat = analyses.tcStatistics(ego_tc, 3, 20);
        
        % put everything into cell arrays (10x10)
        HD_TC{xx,yy} = hd_tc; ALLO_TC{xx,yy} = allo_tc; EGO_TC{xx,yy} = ego_tc;
        HD_ST{xx,yy} = hd_tcStat; ALLO_ST{xx,yy} = allo_tcStat; EGO_ST{xx,yy} = ego_tcStat;
        
    end
end

end

