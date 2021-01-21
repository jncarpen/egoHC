function [model] = modelMe(P, ST, head_direction, params)
%   P:                  position [t x1 y1 x2 y2]
%   ST:                 spiketimes 
%   head_direction:     head direction
%   params:             struct containing initial values for 4 parameters
%   J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET EVERYTHING UP
% parse position data
t = P(:,1); % time (seconds)
fs = mode(diff(t)); % sampling freq
x = P(:,2); 
y = P(:,3);
% x2 = P(:,4);
% y2 = P(:,5);

% spike times (s)
tSpk = ST; 

% remove spike times that are outside the range of tracking times
startTime = t(1); stopTime = t(end);
% histogram edges (binsize close to video frame rate)
edgesT = linspace(startTime,stopTime,numel(t)+1);

% [unsmoothed] spike train
SpkTrn = histcounts(tSpk,edgesT);

tSpk = tSpk(tSpk < stopTime & tSpk > startTime);

% histogram edges (binsize close to video frame rate)
edgesT = linspace(startTime,stopTime,numel(t)+1);

% [unsmoothed] spike train
SpkTrn = histcounts(tSpk,edgesT);

% divide the arena into 100 2D spatial bins
nBins = 10;
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);
yEdges = fliplr(yEdges); % flip y-vector

% get bin centers
for i = 1:length(xEdges)
    if i+1 <= length(xEdges)
        xCenter(i) = ((xEdges(i+1)-xEdges(i))/2)+xEdges(i);
        yCenter(i) = ((yEdges(i+1)-yEdges(i))/2)+yEdges(i);
    end
end

% bin centers (vectors)
count = 1;
for col = 1:length(xCenter)
    for row = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(row), yCenter(col)];
        count = count+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE RATE MAPS (DATA)
% define bins
nBins = 10; % spatial bin size
numBins_HD = 10; % angular bins (36 deg/ea)
binWidth_deg = 360/numBins_HD+1;
angBins = linspace(0,360,numBins_HD+1);

% get angular bin centers
angBinCtrs = [];
for i = 1:length(angBins)
    if i+1 <= length(angBins)
        angBinCtrs(i) = ((angBins(i+1)-angBins(i))/2)+angBins(i);
    end
end


% clear relevant variables
clear loc3d spatial_occupancy rateMap r_xyh R_xyh

% initialize variables
loc3d = zeros(nBins,nBins, 2).* NaN;
spatial_occupancy = zeros(nBins,nBins).*NaN;
rateMap = zeros(nBins,nBins).*NaN;
r_xyh = zeros(nBins,nBins, numBins_HD).*NaN;
R = zeros(nBins,nBins, numBins_HD).*NaN;

% set initial values for the loop
count = 1;
backwardY = nBins:-1:1; forwardX = 1:1:nBins;

% iterate through every 2D spatial bin 
for row = 1:nBins
    for col = 1:nBins
        % what spatial bin are we iterating over?
        yNow = backwardY(row); xNow = forwardX(col);
        
        binNowX(row,col) = xNow;
        binNowY(row,col) = yNow;
        
        % find frames when animal occupied this spatial bin
        indices = find(yNow == binY & xNow == binX);
        timeInBin = length(indices)*fs; % occupancy (s)
        
        % coordinates (in cm) of the location of this particular bin
        loc3d(row,col,:) = [binCenters(count,1), binCenters(count,2)];

        % save occupancy 
        spatial_occupancy(row,col) = timeInBin;

        % spike train for this spatial bin
        spikes_here = SpkTrn(indices); hd_here = head_direction(indices); 
        
        % find head direction (deg) at times of spikes
        spikeInds = find(spikes_here > 0); % since it might be smoothed
        angSpk = hd_here(spikeInds); % angles (in deg)
        % @todo- account for timestamps that have a spikecount > 1

        % make spatial ratemap that does not exclude any bins (pre-threshold)
        % this is for plotting purposes only
        if ~isempty(spikeInds)
            rateMap_inclusive(row,col) = sum(spikes_here./timeInBin);
        else
            rateMap_inclusive(row,col) = 0;
        end
        
        % compute average firing rate for this spatial bin
        r_xy = sum(spikes_here./timeInBin);
        
        % make spatial ratemap
        if timeInBin >= 1 % 1000 ms
            rateMap(row,col) = r_xy; 
        end
        
            % compute angular occupany in each *HD* bin
            [ang_counts, ~, angBinIdx] = histcounts(mod(hd_here, 360), ...
                angBins); 
            
            % probability distribution
            pdx_hd{row,col} = ang_counts./sum(ang_counts);
            
            % calculate angular occupancy (s) for this *spatial bin*
            ang_occupancy_now = ang_counts .* fs;            
            
            for H = 1:length(angBinCtrs)
                % amount of time animal spent in this HD bin (s)
                time_H = ang_occupancy_now(H);
                
                % find indices when animal occupied this HD bin
                idx_H = find(angBinIdx == H);
                
                if time_H >= 0.1 % 100 ms
                    % spiketimes in bin(x,y,H)
                    spk_H = spikes_here(idx_H);
                    
                    % firing rate of cell in bin(x,y,H)
                    r_xyh_now = sum(spk_H./time_H);
                    r_xyh(row,col,H) = r_xyh_now;

                    % FR ratio, R(x,y,H) --> [this is what we compare with model output]
                    if r_xy == 0 && r_xyh_now == 0
                        R_xyh_now = 0;
                        R(row,col,H) = 0;
                    else
                        R_xyh_now = r_xyh_now./r_xy;
                        R(row,col,H) = R_xyh_now;
                    end
                    
                end
            end 
    end
    
    count = count+1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE PARAMETERS & PERFORM OPTIMIZATION

% grab X/Y positions (cm)
X = loc3d(:,:,1); 
Y = loc3d(:,:,2); 

% head direction bins (deg)
H = zeros(nBins, nBins, nBins);
for h = 1:length(angBinCtrs)
    H(:,:,h) = repmat(angBinCtrs(h),10,10);
end

% perform optimization
[output] = fit_jo('cosErr',params,{'g','thetaP','xref','yref'},X,Y,H,R);

% output parameters variable
OP = output.params;

% determine whether the reference point was distant
if OP.xref > 500 || OP.yref > 500
    model.where = 'distant';
else 
    model.where = 'nearby';
end

% correct for negative values of preferred theta(if there are any)
output.params.thetaP = mod(OP.thetaP, 360);

% best fit
[model.pred, model.err] = cosFit(OP, X, Y, H, R);

% get hd tuning curves for each spatial bin (for predicted)
clear pred_values_reshaped
for row=1:nBins
    for col=1:nBins
        pred_values_reshaped{row,col}(:) = model.pred(row,col,:);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE VARIANCE EXPLAINED
%% (1) Variance Explained by Place
var_place = 1 - (var(r_xyh - rateMap, 1, [3 2 1], 'omitnan')...
    ./var(r_xyh, 1, [3 2 1], 'omitnan'));


%% (2) Variance Explained by RH Angle (Model)
var_model = 1 - (var(r_xyh - model.pred, 1, [3 2 1], 'omitnan') ...
    ./var(r_xyh, 1, [3 2 1], 'omitnan'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIND TUNING/MODULATION STRENGTH
% get angular bin centers (rad)
angBinCtrs_rad = deg2rad(angBinCtrs)';

% clear relevant variables
clear mu_rad mu_deg MVL

warning('off','all')
for row = 1:nBins
    for col = 1:nBins
        % grab hd tuning curve in each spatial bin
        tuningCurve_now = reshape(R(row,col,:), 10, 1);
        tuningCurve_now_RH = reshape(model.pred(row,col,:), 10, 1);
        
        % take circular mean (in RADIANS)
        % note: in jercog paper they sum (which would just be mu_rad * 10)
        [mu_rad_uncorrected, ~, ~] = circ_mean(angBinCtrs_rad, ...
            tuningCurve_now);
        mu_rad(row,col) = mod(mu_rad_uncorrected, 2*pi);
        
        [mu_rad_RH_uncorrected, ~, ~] = circ_mean(angBinCtrs_rad, ...
            tuningCurve_now_RH);
        mu_rad_RH(row,col) = mod(mu_rad_RH_uncorrected, 2*pi);
        
        % get circular mean in deg
        mu_deg(row,col) = mod(rad2deg(mu_rad(row,col)), 360);
        mu_deg_RH(row,col) = mod(rad2deg(mu_rad_RH(row,col)), 360);
        
        % mean vector length
        MVL(row,col) = circ_r(angBinCtrs_rad, tuningCurve_now);
        MVL_RH(row,col) = circ_r(angBinCtrs_rad, tuningCurve_now_RH);

        
    end
end
warning('on','all')

% linear average of tuning strengths
MVL(MVL==Inf) = NaN; % MVL(MVL==0) = NaN;
tuningStrength_HD = mean(reshape(MVL, 100,1), 'all', 'omitnan');

MVL_RH(MVL_RH==Inf) = NaN; % MVL_RH(MVL_RH==0) = NaN;
tuningStrength_RH = mean(reshape(MVL_RH, 100,1), 'all', 'omitnan');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE OUTPUTS TO A STRUCT
model.data = R;
model.loc3d = loc3d;
model.datacell = R;
model.bins = angBinCtrs;
model.rateMap = rateMap;
model.rateMapInclusive = rateMap_inclusive;
model.spatbins = binCenters;
model.spatBinNum.X = binNowX;
model.spatBinNum.Y = binNowY;
model.spatial_occ = spatial_occupancy;
model.fval = output.fval;
model.exitflag = output.exitflag;
model.output = output.out;
model.bestParams = output.params;
model.saved = output.saved;
model.varExplained.place = var_place;
model.varExplained.model = var_model;
model.rxyh = r_xyh;
model.modStrength.HD = tuningStrength_HD;
model.modStrength.RH = tuningStrength_RH;
model.modStrength.HD_prefVec = mu_deg;
model.modStrength.RH_prefVec = mu_deg_RH;
model.modStrength.HD_MVL = MVL;
model.modStrength.RH_MVL = MVL_RH;
model.predcell = pred_values_reshaped; % get hd tuning curves for each spatial bin (for predicted)
end



%%
%% option 1: apply a speed threshold
% % speed
% [s, ~] = get_speed(P);
% s = s(:,1); % grab first column
% 
% % speed threshold before we make the spiketrain
% % Get speed at time of spike and put into vector SpikeSpeed
% SpikeSpeed = interp1 (t, s, ST); %in cm/s
% 
% % Set threshold
% thr_d= 4; % this is the threshold set in jercog et al. (diff for dMan)       
% thr_u= 100;
% 
% % Apply threshold 
% a=find(SpikeSpeed>thr_d);
% b=find(SpikeSpeed<thr_u);
% 
% % make position samples NaN
% x(find(s<thr_d))=NaN; x(find(s>thr_u))=NaN;
% y(find(s<thr_d))=NaN; y(find(s>thr_u))=NaN;
% x2(find(s<thr_d))=NaN; x2(find(s>thr_u))=NaN;
% y2(find(s<thr_d))=NaN; y2(find(s>thr_u))=NaN;
% position = [t, x, y, x2, y2];
% 
% % Combined threshold 
% c=intersect(a,b);
% 
% % Vector with filtered spikes - based on indexing from c
% SpikeSpeed_fil=ST(c);
% tSpk = SpikeSpeed_fil; % spike times
% 
% % MAKE SPIKE TRAIN (bin the spikes)- this is speed-thresholded
% startTime = t(1);
% stopTime = t(end);
% 
% % remove spike times that are outside the range of tracking times
% tSpk = tSpk(tSpk < stopTime & tSpk > startTime);
% 
% edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate
% 
% binnedSpikes = histcounts(tSpk,edgesT);
% 
% sigma = 2; % smoothing factor
% SpkTrn = imgaussfilt(binnedSpikes, sigma, 'Padding', 'replicate'); % smooth spiketrain
% 
% % get head direction values
% % head_direction = get_hd(position);
% 
% % make sure that head direction ranges from 0-360 deg
% if nanmax(head_direction) < 350
%     head_direction = rad2deg(head_direction);
% else
%     head_direction = head_direction;
% end

