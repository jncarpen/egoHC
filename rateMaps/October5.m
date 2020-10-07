% Messy script for lab meeting presentation
% October 5, 2020

% WHOS DATA ARE WE LOOKING AT?:

%% (1) for Jan Sigurd's data:
% define variables for the session you're interested in
sessNum = 27; unitNum = 4;
position = pos_cm{1,sessNum}; head_direction = hd{1,sessNum};
[spkTrainNow_thresh] = spkTrain_thresh(pos_cm, SpikeTimes_thresh, SpikeTrain); % get speed-thresholded spiketrain
SpkTrn = spkTrainNow_thresh{1,sessNum}{1,unitNum}';


%% (2) for Sebastian's data:
% correct position data (cm) 
boxSize = 80; % in cm (according to Seb)
xLen = nanmax(positions(:,2))-nanmin(positions(:,2));
xLen2 = nanmax(positions(:,4))-nanmin(positions(:,4));
yLen = nanmax(positions(:,3))-nanmin(positions(:,3));
yLen2 = nanmax(positions(:,5))-nanmin(positions(:,5));
conFac_x = boxSize/xLen; conFac_y = boxSize/yLen;
conFac_x2 = boxSize/xLen2; conFac_y2 = boxSize/yLen2;
t = positions(:,1); sampleRate = mode(diff(t));
x = (positions(:,2)-nanmin(positions(:,2)))*conFac_x; y = (positions(:,3)-nanmin(positions(:,3)))*conFac_y;
x2 = (positions(:,4)-nanmin(positions(:,4)))*conFac_x2; y2 = (positions(:,5)-nanmin(positions(:,5)))*conFac_y2;
position = [t, x, y, x2, y2]; % corrected position

% get head_direction values
head_direction = rem(atan2d(position(:,5)-position(:,3),position(:,4)-position(:,2)) + 180, 360);

% speed threshold spikes
t = position(:,1);
x_smooth=medfilt1(x);y_smooth=medfilt1(y);
speed_OVC = zeros(length(t), 1);
for i = 2:numel(x_smooth)-1
    speed_OVC(i) = sqrt((x_smooth(i+1) - x_smooth(i-1))^2 + (y_smooth(i+1) - y_smooth(i-1))^2) / (t(i+1) - t(i-1));
end

% make a spike train
stopTime = nanmax(t); startTime = nanmin(t);
S = cellTS(cellTS < stopTime & cellTS > startTime); % remove times outside of recording
edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate
binnedSpikes = histcounts(S,edgesT); % bin those spikes baby
speed_idx = find(speed_OVC<5); % find indices when animal was moving slow
binnedSpikes(speed_idx)=0; % get rid of spikes when animal was moving slow
binnedSpikes = imgaussfilt(binnedSpikes, 2, 'Padding', 'replicate'); % smooth ST
SpkTrn = binnedSpikes;


%% (3) Sebastian's data (simulated EBC)
% find egobearing for each timepoint
rlX=42; rlY=59; % object moved
midX=(x+x2)/2; midY=(y+y2)/2;
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
egoAng = alloAng - head_direction;

% define angles of interest
angle_of_interest = 60; plus_minus = 20;
min_angle = angle_of_interest - plus_minus;
max_angle = angle_of_interest + plus_minus;
logical = egoAng>min_angle & egoAng<max_angle; % find indices where egoAngle is within range
idx = find(logical==1); %logical = alloAng>min_angle & alloAng<max_angle;
throw_away = .2; % percentage to remove
sz = floor(length(idx)-length(idx)*throw_away); % size to keep
len_background = sz*.3; % number of background spikes to add
randIdx = datasample(idx, sz, 'Replace', false);
foreground_spikes = t(randIdx);
simTS = sort(foreground_spikes, 'ascend'); % simulated timestamps

% speed threshold spikes
t = position(:,1);
x_smooth=medfilt1(x);y_smooth=medfilt1(y);
speed_OVC = zeros(length(t), 1);
for i = 2:numel(x_smooth)-1
    speed_OVC(i) = sqrt((x_smooth(i+1) - x_smooth(i-1))^2 + (y_smooth(i+1) - y_smooth(i-1))^2) / (t(i+1) - t(i-1));
end

% make a spike train
stopTime = nanmax(t); startTime = nanmin(t);
S = simTS(simTS < stopTime & simTS > startTime); % remove times outside of recording
edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate
binnedSpikes = histcounts(S,edgesT); % bin those spikes baby
speed_idx = find(speed_OVC<5); % find indices when animal was moving slow
binnedSpikes(speed_idx)=0; % get rid of spikes when animal was moving slow
binnedSpikes = imgaussfilt(binnedSpikes, 2, 'Padding', 'replicate'); % smooth ST
SpkTrn = binnedSpikes;


%% run the analysis
% get a vector of all the reference points you're interested in
nBins = 10; % increase this # when ready to run
[refVec, outside_circ] = generate_reference_pnts(position, "True", nBins); % the number of reference points you get is nBins^2

% loop through vector of reference points (refVec)
for refPt = 1:nBins^2
    % grab reference point for this iteration
    refLoc = refVec(refPt,:);
    
    % compute tuning curves + stats for this *reference point*
    [HD_TC, ALLO_TC, EGO_TC, HD_ST, ALLO_ST, EGO_ST] = TC_stats_2DBins(position, head_direction, SpkTrn, refLoc);
    
    % get scores for this reference points
    [HD_mean_stats, HD_sum_stats] = score_tuning_curve(HD_ST);
    [ALLO_mean_stats, ALLO_sum_stats] = score_tuning_curve(ALLO_ST);
    [EGO_mean_stats, EGO_sum_stats] = score_tuning_curve(EGO_ST);
    
    % lets test MVL
    ego_test_MVL(refPt) = EGO_sum_stats.MVL;
    allo_test_MVL(refPt) = ALLO_sum_stats.MVL;
    hd_test_MVL(refPt) = HD_sum_stats.MVL;
    
    ego_test_MD(refPt) = EGO_mean_stats.mean_direction;
    allo_test_MD(refPt) = ALLO_mean_stats.mean_direction;
    hd_test_MD(refPt) = HD_mean_stats.mean_direction;
    
    ego_test_pr(refPt) = EGO_sum_stats.peak_rate;
    allo_test_pr(refPt) = ALLO_sum_stats.peak_rate;
    hd_test_pr(refPt) = HD_sum_stats.peak_rate;
    
    
end

test2 = ego_test_MVL;
% reshape the matrix
start = 1; stop = 100;
for ii = 1:100
%     meanMat(:,ii) = flip(test(start:stop))';
    sumMat(:,ii)= flip(test2(start:stop))';
    start = start + nBins; stop = stop + nBins;
end

% plot results (sum mat)
figure
set(gcf,'color','w');
imagesc(sumMat)
title("S27U4:HD")
pbaspect([1 1 1])
% colormap(hsv)
% caxis([0 360])
colorbar
% set(gca,'YDir','normal')
box off

% % plot results (mean mat)
% figure
% set(gcf,'color','w');
% imagesc(meanMat)
% title("EBC-Object trial data")
% pbaspect([1 1 1])
% colorbar
% % set(gca,'YDir','normal')
% box off

%% plot tuning curves for 



