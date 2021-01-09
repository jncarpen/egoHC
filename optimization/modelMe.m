function [model] = modelMe(P, ST, head_direction, params)
%   P:                  position [t x1 y1 x2 y2]
%   ST:                 spiketimes 
%   head_direction:     head direction
%   params:             struct containing initial values for 4 parameters
%   J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get information about cell

% parse position data
t = P(:,1); % time (seconds)
fs = mode(diff(t)); % sampling freq
x = P(:,2); 
y = P(:,3);
x2 = P(:,4);
y2 = P(:,5);

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

%% option 2: make a spiketrain without speed thresholding
tSpk = ST; % spike times
% position = P; % rename position

% MAKE SPIKE TRAIN (bin the spikes)- this is speed-thresholded
startTime = t(1);
stopTime = t(end);

% remove spike times that are outside the range of tracking times
tSpk = tSpk(tSpk < stopTime & tSpk > startTime);

edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate

% bin the spikes
binnedSpikes = histcounts(tSpk,edgesT);

% smooth the spike train (by some factor, sigma)
sigma = 2; % smoothing factor
SpkTrn = imgaussfilt(binnedSpikes, sigma, 'Padding', 'replicate'); % smooth spiketrain


%% bin the spatial arena
% divide the arena into 100 2D spatial bins
nBins = 10;
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);
yEdges = fliplr(yEdges); % flip y-vector

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

count = 1;
for col = 1:length(xCenter)
    for row = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(row), yCenter(col)];
%         plot(xCenter(row), yCenter(col), '.'); hold on;
        count = count+1;
    end
end

%% compute values for each bin

% define spatial bins
numBins_HD = 10; % 24 degree bins
binWidth_deg = 360/numBins_HD+1;
angBins = linspace(0,360,numBins_HD+1);

% get angular bin centers
angBinCtrs = [];
for i = 1:length(angBins)
    if i+1 <= length(angBins)
        angBinCtrs(i) = ((angBins(i+1)-angBins(i))/2)+angBins(i);
    end
end

% clear old variables
clear rateMap rateMap_HD R_ratio R_ratio_summed pdx_hd r_xyh_mat ...
    ang_occupancy

% initialize variables
spatial_occupancy=zeros(nBins,nBins);
RR = zeros(nBins, nBins, numBins_HD);
R_ratio = cell(nBins,nBins);
angOccArray = cell(nBins,nBins);

nBins = 10;
count = 1;
didNotPass = 0;
backwardY = nBins:-1:1; forwardX = 1:1:nBins;


% iterate through every 2D spatial bin 
for row = 1:nBins
    for col = 1:nBins
        yNow = backwardY(row);
        xNow = forwardX(col);
        indices = find(yNow == binY & xNow == binX);
        timeInBin = length(indices)*fs; % occupancy (s)
        
        % x and y now in matrix form 
        xNow_(row,col) = xNow;
        yNow_(row,col) = yNow;
        
        % save the spatial bin numbers for each count
        spatbinNum(count,1) = xNow;
        spatbinNum(count,2) = yNow;
        
        % coordinates (in cm) of the location of this particular bin
        location_now{row,col} = [binCenters(count,1), binCenters(count,2)]; % array form
        loc3d(row,col,:) = [binCenters(count,1), binCenters(count,2)]; % matrix form

        % save occupancy 
        spatial_occupancy(row,col) = timeInBin;

        % spike train for this *spatial* bin
        spikes_here = SpkTrn(indices); hd_here = head_direction(indices); 
        
%         disp(strcat('bin ', sprintf('%.f', count), 'has ', sprintf('%.2f', sum(spikes_here)), ' spikes...'))
        
        % find head direction at times of spikes
        spikeInds = find(spikes_here > 0); % since its smoothed
        angSpk = hd_here(spikeInds); % angles (in deg)

        % make ratemap that does not exclude any bins (pre-threshold)
        if ~isempty(spikeInds)
            rateMap_inclusive(row,col) = sum(spikes_here./timeInBin);
        else
            rateMap_inclusive(row,col) = 0;
        end
        
        % set a bin-by-bin threshold
        if ~isempty(spikeInds) && length(spikeInds) > 10 && timeInBin > 0.5
        % if ~isempty(spikeInds) && length(spikeInds) > 5 && timeInBin > 0.5 && all(ang_occupancy > .05)
        % are at least 50 degrees covered?
        
            % compute angular occupany in each *HD* bin
            [ang_counts, ~, angBinIdx] = histcounts(mod(hd_here, 360), angBins); % circular histogram
            pdx_hd{row,col} = ang_counts./sum(ang_counts); % probability distribution
            
            % check to make sure the pdx sums to 1
%             if not(sum(pdx_hd{row,col}, 'omitnan') == 1)
%                 disp(strcat('pdx_hd', '(', sprintf('%.f', row), ',' ,sprintf('%.f', col), ') does not sum to 1.'));
%             end
            
            % calculate angular occupancy for this *spatial bin*
            ang_occupancy(row,col, :) = ang_counts .* fs; % how many seconds the animal was in each angular bin
            
            
            % compute normal firing rate map [r(x,y)] for this 2D spatial bin
            r_xy = sum(spikes_here./timeInBin); % firing rate in this spatial bin
            rateMap(row,col) = r_xy; % save for later
            
            for H = 1:length(angBinCtrs)
                % find indices of timestamps when animal occupied both
                % spatial_bin(row,col) and HD_bin(H) --> (bin(x,y,H))
                idx_H = find(H == angBinIdx);
                
                % amount of time animal spend in bin(x,y,H)
                time_H = length(idx_H)*fs;
                
                % spiketimes in bin(x,y,H) 
                spk_H = spikes_here(idx_H);
                
                % save angular [spike] occupancy (in counts)
                spk_H(spk_H>0) = 1;
                angOccArray{row,col}(H,1) = sum(spk_H, 'omitnan');
                
                % firing rate of cell in bin(x,y,H)
                r_xyh = sum(spk_H./time_H);
                rateMap_HD{row,col}(H,1) = r_xyh; % save in cell array
                r_xyh_mat(row,col,H) = r_xyh; % save in matrix
                
                % FR ratio in bin(x,y,H), R(x,y,H) --> [this is what we compare with model output]
                R_xyh = r_xyh./r_xy; 
                R_ratio{row,col}(H,1) = R_xyh; % save in cell array 
                RR(row,col,H) = R_xyh; % save in matrix
            
            end
            
            % summed R_ratio, @todo: comment this better later... lol
            R_ratio_summed(row,col) = nansum(R_ratio{row,col});

        else
           % disp(strcat('bin ', sprintf('%.f', count), ' did not pass criteria...'))
           didNotPass = didNotPass + 1; % add one to the 'didNotPass' count
           rateMap(row,col) = NaN; % throw a NaN into the rateMap
           r_xyh_mat(row,col,:) = NaN; % throw a NaN into the rateMap that is conditioned on HD
           angOccArray{row,col} = [];
           ang_occupancy(row,col, :) = zeros(1,10)*NaN;
        end
        
        count = count+1;
    end
end

% display number of spatial bins that didn't pass the criteria
% disp(didNotPass)


%% find HD tuning strength for this neuron (real data)
% get angular bin centers in radians
angBinCtrs_rad = deg2rad(angBinCtrs)';

clear mu_rad mu_deg MVL
warning('off','all')
for row = 1:nBins
    for col = 1:nBins
        % grab hd tuning curve in each spatial bin
        tuningCurve_now = reshape(RR(row,col,:), 10, 1);
        
        % get angular occupancy
        angOcc_now = angOccArray{row,col};
        
        % take circular mean (in RADIANS)
        % note: in jercog paper they sum (which would just be mu_rad * 10)
        [mu_rad_uncorrected, ~, ~] = circ_mean(angBinCtrs_rad, tuningCurve_now);
        mu_rad(row,col) = mod(mu_rad_uncorrected, 2*pi);
        
        % get circular mean in deg
        mu_deg(row,col) = mod(rad2deg(mu_rad(row,col)), 360);
        
        
        if ~isempty(angOcc_now)
            % find mean resultant length
            MVL(row,col) = circ_r(angBinCtrs_rad, angOcc_now);
        else
            MVL(row,col) = 0;
        end
        
    end
end
warning('on','all')

% take mean of MVLs to get head direction tuning strength
% for this unit
MVL(MVL==Inf) = NaN; % get rid of infinity values
MVL(MVL==0) = NaN;
tuningStrength_HD = mean(reshape(MVL, 100,1), 'all', 'omitnan');


%% fit the model

% get center of the box;
% [boxCtrX,boxCtrY] = getBoxCenter(P);

% set initial values for parameters to be fit
% clear p
% p.g = .25;
% p.thetaP = 0;
% p.xref = boxCtrX;
% p.yref = boxCtrY;

clear p
p.g = params.g;
p.thetaP = params.thetaP;
p.xref = params.xref;
p.yref = params.yref;

% set values for stable parameters
X = loc3d(:,:,1); 
Y = loc3d(:,:,2); 

H = zeros(nBins, nBins, nBins);
for h = 1:length(angBinCtrs)
    H(:,:,h) = repmat(angBinCtrs(h),10,10);
end

% make a copy of RR (save as R)
R = RR;

% have the model take a first guess given the parameters we've fed in
% [firstGuess_pred, firstGuess_err] = cosFit(p,X,Y,H,R);

% fit the model (make sure this is the modified 'fit' function (not the orig))
clear output
[output] = fit_jo('cosErr',p,{'g','thetaP','xref','yref'},X,Y,H,R);

% reassign variable name (debugging)
OP = output.params;

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


%% compute variance explained (as in Jercog et al.)
% (1) variance explained by place tuning (?)
num = r_xyh_mat - rateMap;
num = var(num, 1, [3 2 1], 'omitnan');
den = var(r_xyh_mat, 1, [3 2 1], 'omitnan');
var_place = 1 - (num/den);

% (2) variance explained overall (by the model)
numMod = r_xyh_mat - model.pred;
numMod = var(numMod, 1, [3 2 1], 'omitnan');
denMod = var(r_xyh_mat, 1, [3 2 1], 'omitnan');
var_model = 1 - (numMod/denMod);


%% find RH tuning strength for this neuron (model data)
% get angular bin centers in radians
angBinCtrs_rad = deg2rad(angBinCtrs)';

clear mu_rad_RH mu_deg_RH MVL_RH
warning('off','all')
for row = 1:nBins
    for col = 1:nBins
        % grab hd tuning curve in each spatial bin
        tuningCurve_now = reshape(model.pred(row,col,:), 10, 1);
%         tuningCurve_now(tuningCurve_now<0)=0;
        
        % take circular mean (in RADIANS)
        % note: in jercog paper they sum (which would just be mu_rad * 10)
        [mu_rad_RH_uncorrected, ~, ~] = circ_mean(angBinCtrs_rad, tuningCurve_now);
        mu_rad_RH(row,col) = mod(mu_rad_RH_uncorrected, 2*pi);
        
        
        % get circular mean in deg
        mu_deg_RH(row,col) = mod(rad2deg(mu_rad_RH(row,col)), 360);
        
        % find mean resultant length
        MVL_RH(row,col) = circ_r(angBinCtrs_rad, tuningCurve_now);
        
    end
end
warning('on','all')

% take mean of MVLs to get head direction tuning strength
% for this unit
MVL_RH(MVL_RH==0) = NaN;
MVL_RH(MVL_RH==Inf) = NaN; % get rid of infinity values
tuningStrength_RH = mean(reshape(MVL_RH, 100,1), 'all', 'omitnan');


%% save outputs in a struct
model.data = R;
model.loc3d = loc3d;
model.datacell = R_ratio;
model.bins = angBinCtrs;
model.rateMap = rateMap;
model.rateMapInclusive = rateMap_inclusive;
model.spatbins = binCenters;
model.spatbinsnum = spatbinNum;
model.spatial_occ = spatial_occupancy;
model.fval = output.fval;
model.exitflag = output.exitflag;
model.output = output.out;
model.bestParams = output.params;
model.saved = output.saved;
model.varExplained.place = var_place;
model.varExplained.model = var_model;
model.rxyh = r_xyh_mat;
model.modStrength.HD = tuningStrength_HD;
model.modStrength.RH = tuningStrength_RH;
model.modStrength.HD_prefVec = mu_deg;
model.modStrength.RH_prefVec = mu_deg_RH;
model.modStrength.HD_MVL = MVL;
model.modStrength.RH_MVL = MVL_RH;
model.predcell = pred_values_reshaped; % % get hd tuning curves for each spatial bin (for predicted)
end

