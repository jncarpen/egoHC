% Define variables and discreteize to feed into a GLM 
% (Originally written for the OVC project- still needs to be adapted)
% and simplified...

% Jordan Carpenter 
% Last Modified: May 6, 2020

% initialize empty cell array to store speed values for all conditions
speed = cell(1,3);
PosX = cell(1,3);
PosY = cell(1,3);
ObjDistance = cell(1,3);
HeadDirection = cell(1,3);
smoothSpikes = cell(1,3);
sqrtSpikes = cell(1,3);
TrackingTimes = cell(1,3);
hdDiscretize = cell(1,3);
speedDiscretize = cell(1,3);
disDiscretize = cell(1,3);
posXDiscretize = cell(1, 3);
posYDiscretize = cell(1,3);
spiketrain = cell(1,3);
% thetaPhase = cell(1,3);
dataset = cell(1,7);

% define object coordinates 
ObjCoord_2 = coordinates(1,1:2);% coord for condition #2
ObjCoord_3 = coordinates(2,1:2);% coord for condition #3


% loop through all three conditions
for condition = 1:3
%%  make poisson spikes look gaussian

% pull out data for condition
PosC = Pos{1,condition};
EEGC = nonzeros(EEG{1, condition}{1,1});
spikeIndsC = spikeInds{1, condition}; % scale to match time vector (seconds)

% define time vector for condition
trackingtimes = PosC(:,1); % in seconds

% save trackingtimes for each condition
TrackingTimes{1, condition} = trackingtimes;

% compute first and last timestamp
startTime = trackingtimes(1);
stopTime = trackingtimes(end);

% define sample rate
sampleRate = mode(diff(trackingtimes));

% intiialize empty cell array to store smoothed spike info for each cell
smoothCell = cell(1, length(spikeIndsC));
pSpikes = cell(1, length(spikeIndsC));
% sqrtCell = cell(1, length(spikeIndsC));

for neuron = 1:length(spikeIndsC) % loop through all cells for the condition
    
    % get spikes for current cell
    spikeIdxPos = spikeIndsC{1,neuron}(spikeIndsC{1,neuron} ~= 101); % remove excluded spikes
    spikeTimes = Pos{1,condition}(spikeIdxPos,1);
    
    % remove spike times that are outside the range of tracking times
    spikeTimes = spikeTimes(spikeTimes < stopTime & spikeTimes > startTime);
    
    % bin spike data (by time)
    edgesT = linspace(startTime,stopTime,numel(trackingtimes)+1); % binsize is close to video frame rate
    
    % histcount of spikes per time bin
    binnedSpikes = histcounts(spikeTimes,edgesT);
    
    % compute bin size for gaussian window
    %binSz = round_odd(std(trackingtimes)*3);
    binSz = 5;
    
    % smooth time-binned data w/ a gaussian filter (what Ben said)
    binnedSpikes_smooth = smoothdata(binnedSpikes, 'gaussian', binSz);
    binnedSpikes_smooth = sqrt(binnedSpikes_smooth)';
    
    % store smoothed values for each neuron in a cell array
    smoothCell{1, neuron} = binnedSpikes_smooth;
    pSpikes{1, neuron} = binnedSpikes';
    %sqrtCell{1, neuron} = binnedSpikes;
    
end
    
    % store smoothed values for each condition in a cell array
    smoothSpikes{1, condition} = smoothCell;
    spiketrain{1, condition} = pSpikes;
    % sqrtSpikes{1, condition} = sqrtCell;
    
%% compute speed

% define posx, posy, posx2, posy2
posx = inpaint_nans(PosC(:,2),3); % interpolate nan values w/ a del^4 interpolating operator
posy = inpaint_nans(PosC(:,3),3);
posx2 = inpaint_nans(PosC(:,4),3);
posy2 = inpaint_nans(PosC(:,5),3);

% define x-coordinates for LED1
%posxMean = smoothdata(nanmean(PosC(:,[2 4]), 2),'gaussian', 7); % returns the mean, or uses one of the x values if the other is NaN
posxMean = smoothdata(posx, 'gaussian', binSz); 

% define y-coordinates for LED1
%posyMean = smoothdata(nanmean(PosC(:,[3 5]), 2),'gaussian', 7);
posyMean = smoothdata(posy, 'gaussian', binSz); 

% compute distance between points for x and y coordinates
% this is where I get zero values (crossover points?)//
% the LEDs are switching at certain points *
velx = diff([posxMean(1); posxMean]); vely = diff([posyMean(1); posyMean]); % add the extra just to make the vectors the same size

% output unsmoothed speed
Spd = (sqrt(velx.^2+vely.^2).*sampleRate)*1000;
Spd = medfilt1(Spd, 7);
% Spd = smoothdata(Spd, 'gaussian', 7);

% save relevant variables in cell arrays
speed{1, condition} = Spd;
PosX{1, condition} = posxMean;
PosY{1, condition} = posyMean;

%% discretize position

% bin position data
posBins = 20;
posEdges = linspace(0, 80, posBins+1);

%%%%%%%%%%%%%%%%%%%%%%%

% For X-position:

% compute position occupancy & pos indices
[posXOccupancy,~,posXInds] = histcounts(posxMean, posEdges);

% initialize zero matrix
posX_dis = zeros(length(posXInds), posBins);

for T=1:length(posXInds) % for every time point
    currentBin = posXInds(T); % pulls out value for HD bin at time_i
    posX_dis(T, currentBin+1) = 1; % sets value to 1
end

posXDiscretize{1, condition} = posX_dis;

%%%%%%%%%%%%%%%%%%%%%%%

% For Y-position:

% compute position occupancy & pos indices
[posYOccupancy,~,posYInds] = histcounts(posyMean, posEdges);

% initialize zero matrix
posY_dis = zeros(length(posYInds), posBins);

for T=1:length(posYInds) % for every time point
    currentBin = posYInds(T); % pulls out value for HD bin at time_i
    posY_dis(T, currentBin+1) = 1; % sets value to 1
end

posYDiscretize{1, condition} = posY_dis;


%% discretize speed

% bin speed data
speedBins = 20;
speedEdges = linspace(0,25,speedBins+1);

% compute speed occupancy & speed indices
[speedOccupancy,~,speedInds] = histcounts(Spd, speedEdges);

% initialize zero matrix
speed_dis = zeros(length(speedInds), speedBins);

for T=1:length(speedInds) % for every time point
    currentBin = speedInds(T); % pulls out value for HD bin at time_i
    speed_dis(T, currentBin+1) = 1; % sets value to 1
end

speedDiscretize{1, condition} = speed_dis;

%% compute head direction

% compute head direction
direction = atan2(posy2-posy,posx2-posx)+pi/2;
direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers

% bin HD data 
n_bins_angle = 20;
edgesHD = linspace(0,2*pi,n_bins_angle+1);

% compute occupancy & angle indices
[occupancy,~,angle_inds] = histcounts(direction,edgesHD);


% save HD for each condition in a cell array
HeadDirection{1, condition} = direction;

%% discretize HD
% initialize zero matrix
HD_dis = zeros(length(angle_inds), n_bins_angle);

for T=1:length(angle_inds) % for every time point
    currentBin = angle_inds(T); % pulls out value for HD bin at time_i
    HD_dis(T, currentBin) = 1; % sets value to 1
end

hdDiscretize{1, condition} = HD_dis;
%% Visualization: plot occ and map
% for i=1:100
% for j= 1:100
% f = find(abs(posxMean-i)<5 & abs(posyMean-j)<5);
% occ(i,j) = length(f);
% map(i,j) = nanmedian(Spd(f));
% end
% end
% 
% figure 
% imagesc(occ)
% title("Occupancy")
% % xlabel("x")
% % ylabel("y")
% 
% figure
% imagesc(map)
% caxis([.0025 .015]);
% title("Speed Map")
% % xlabel("x")
% % % ylabel("y")

%% compute distance from object

if condition == 1 % compute distance from object in condition 2
     D = [];
    for cordInd = 1:length(posxMean) % loop through all time positions in posx and posy
        % use distance formula to calculate distance from LED to object
        DisValue = sqrt((ObjCoord_2(1,1)- posxMean(cordInd))^2+ (ObjCoord_2(1,2)- posyMean(cordInd))^2);
        D = [D; DisValue];
    end
end

if condition == 2 % for initial object condition
    D = [];
    for cordInd = 1:length(posxMean) % loop through all time positions in posx and posy
        % use distance formula to calculate distance from LED to object
        DisValue = sqrt((ObjCoord_2(1,1)- posxMean(cordInd))^2+ (ObjCoord_2(1,2)- posyMean(cordInd))^2);
        D = [D; DisValue];
    end
end

if condition == 3
    D = [];
    for cordInd = 1:length(posxMean)
        % use distance formula to calculate distance from LED to object
        DisValue = sqrt((ObjCoord_2(1,1)- posxMean(cordInd))^2+ (ObjCoord_2(1,2)- posyMean(cordInd))^2);
        D = [D; DisValue];
    end
end

ObjDistance{1, condition} = D;

%% discretize distance from object

% bin dist from object data 
disBin = 20;
disEdges = linspace(0,50,disBin+1);

% compute occupancy & angle indices
[disOccupancy,~,disInds] = histcounts(D, disBin);

Dis_dis = zeros(length(disInds), disBin);

for T=1:length(disInds) % for every time point
    currentBin = disInds(T); % pulls out value for HD bin at time_i
    Dis_dis(T, currentBin) = 1; % sets value to 1
end

disDiscretize{1, condition} = Dis_dis;

%% compute instantaneous phase
% deal with zeroes at the end of LFP signal
eegSR = 250; % described in EEG file
[b a] = butter(3,[4/eegSR/2 12/eegSR/2],'bandpass');
filtered = filtfilt(b,a,EEGC);
phase = angle(hilbert(filtered)); % phase angle?

% hilb_eeg = hilbert(EEGC); % compute hilbert transform
% phase = atan2(imag(hilb_eeg),real(hilb_eeg)); %inverse tangent (-pi to pi)
% ind = phase <0; phase(ind) = phase(ind)+2*pi; % from 0 to 2*pi
%
end


%% store cell arrays in master cell array
%dataset{1,1} = "smoothSpikes";
dataset{1,1} = smoothSpikes;
%dataset{1,2} = "PosX";
dataset{1,2} = PosX;
%dataset{1,3} = "PosY";
dataset{1,3} = PosY;
%dataset{1,4} = "HeadDirection";
dataset{1,4} = HeadDirection;
%dataset{1,5} = "speed";
dataset{1,5} = speed;
%dataset{1,6} = "ObjDistance";
dataset{1,6} = ObjDistance;
dataset{1,7} = TrackingTimes;


%% Function: Round Odd
function S = round_odd(S)
% round to nearest odd integer.
idx = mod(S,2)<1;
S = floor(S);
S(idx) = S(idx)+1;
end

%% Function: Interpolate NaNs
function B=inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%         
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%         
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this
%         method is disabled for vector inputs.
%         
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%         
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
  method = 0;
elseif ~ismember(method,0:5)
  error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
 case 0
  % The same as method == 1, except only work on those
  % elements which are NaN, or at least touch a NaN.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really a 1-d case
    work_list = nan_list(:,1);
    work_list = unique([work_list;work_list - 1;work_list + 1]);
    work_list(work_list <= 1) = [];
    work_list(work_list >= nm) = [];
    nw = numel(work_list);
    
    u = (1:nw)';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
      repmat([1 -2 1],nw,1),nw,nm);
  else
    % a 2-d case
    
    % horizontal and vertical neighbors only
    talks_to = [-1 0;0 -1;1 0;0 1];
    neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
    
    % list of all nodes we have identified
    all_list=[nan_list;neighbors_list];
    
    % generate sparse array with second partials on row
    % variable for each element in either list, but only
    % for those nodes which have a row index > 1 or < n
    L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    else
      fda=spalloc(n*m,n*m,size(all_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    end
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 1
  % least squares approach with del^2. Build system
  % for every array element as an unknown, and then
  % eliminate those which are knowns.

  % Build sparse matrix approximating del^2 for
  % every element in A.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % a 1-d case
    u = (1:(nm-2))';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
      repmat([1 -2 1],nm-2,1),nm-2,nm);
  else
    % a 2-d case
    
    % Compute finite difference for second partials
    % on row variable first
    [i,j]=ndgrid(2:(n-1),1:m);
    ind=i(:)+(j(:)-1)*n;
    np=(n-2)*m;
    fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      repmat([1 -2 1],np,1),n*m,n*m);
    
    % now second partials on column variable
    [i,j]=ndgrid(1:n,2:(m-1));
    ind=i(:)+(j(:)-1)*n;
    np=n*(m-2);
    fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
      repmat([1 -2 1],np,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 2
  % Direct solve for del^2 BVP across holes

  % generate sparse array with second partials on row
  % variable for each nan element, only for those nodes
  % which have a row index > 1 or < n
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really just a 1-d case
    error('Method 2 has problems for vector input. Please use another method.')
    
  else
    % a 2-d case
    L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    else
      fda=spalloc(n*m,n*m,size(nan_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    end
    
    % fix boundary conditions at extreme corners
    % of the array in case there were nans there
    if ismember(1,nan_list(:,1))
      fda(1,[1 2 n+1])=[-2 1 1];
    end
    if ismember(n,nan_list(:,1))
      fda(n,[n, n-1,n+n])=[-2 1 1];
    end
    if ismember(nm-n+1,nan_list(:,1))
      fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
    end
    if ismember(nm,nan_list(:,1))
      fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
    end
    
    % eliminate knowns
    rhs=-fda(:,known_list)*A(known_list);
    
    % and solve...
    B=A;
    k=nan_list(:,1);
    B(k)=fda(k,k)\rhs(k);
    
  end
  
 case 3
  % The same as method == 0, except uses del^4 as the
  % interpolating operator.
  
  % del^4 template of neighbors
  talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
      0 1;0 2;1 -1;1 0;1 1;2 0];
  neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
  
  % list of all nodes we have identified
  all_list=[nan_list;neighbors_list];
  
  % generate sparse array with del^4, but only
  % for those nodes which have a row & column index
  % >= 3 or <= n-2
  L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
  nl=length(L);
  if nl>0
    % do the entire template at once
    fda=sparse(repmat(all_list(L,1),1,13), ...
        repmat(all_list(L,1),1,13) + ...
        repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
        repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
  else
    fda=spalloc(n*m,n*m,size(all_list,1)*5);
  end
  
  % on the boundaries, reduce the order around the edges
  L = find((((all_list(:,2) == 2) | ...
             (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
           (((all_list(:,3) == 2) | ...
             (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,5), ...
      repmat(all_list(L,1),1,5) + ...
        repmat([-n,-1,0,+1,n],nl,1), ...
      repmat([1 1 -4 1 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,2) == 1) | ...
             (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-n,0,n],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,3) == 1) | ...
             (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-1,0,1],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 4
  % Spring analogy
  % interpolating operator.
  
  % list of all springs between a node and a horizontal
  % or vertical neighbor
  hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
  hv_springs=[];
  for i=1:4
    hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
    k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
    hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
  end

  % delete replicate springs
  hv_springs=unique(sort(hv_springs,2),'rows');
  
  % build sparse matrix of connections, springs
  % connecting diagonal neighbors are weaker than
  % the horizontal and vertical springs
  nhv=size(hv_springs,1);
  springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
     repmat([1 -1],nhv,1),nhv,nm);
  
  % eliminate knowns
  rhs=-springs(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
  
 case 5
  % Average of 8 nearest neighbors
  
  % generate sparse array to average 8 nearest neighbors
  % for each nan element, be careful around edges
  fda=spalloc(n*m,n*m,size(nan_list,1)*9);
  
  % -1,-1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,-1
  L = find(nan_list(:,3) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,-1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,0
  L = find(nan_list(:,2) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,0
  L = find(nan_list(:,2) < n);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,+1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,+1
  L = find(nan_list(:,3) < m);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,+1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  k=nan_list(:,1);
  B(k)=fda(k,k)\rhs(k);
  
end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);


% ====================================================
%      end of main function
% ====================================================
% ====================================================
%      begin subfunctions
% ====================================================
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%      
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
% 
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
  % use the definition of a neighbor in talks_to
  nan_count=size(nan_list,1);
  talk_count=size(talks_to,1);
  
  nn=zeros(nan_count*talk_count,2);
  j=[1,nan_count];
  for i=1:talk_count
    nn(j(1):j(2),:)=nan_list(:,2:3) + ...
        repmat(talks_to(i,:),nan_count,1);
    j=j+nan_count;
  end
  
  % drop those nodes which fall outside the bounds of the
  % original array
  L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
  nn(L,:)=[];
  
  % form the same format 3 column array as nan_list
  neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
  
  % delete replicates in the neighbors list
  neighbors_list=unique(neighbors_list,'rows');
  
  % and delete those which are also in the list of NaNs.
  neighbors_list=setdiff(neighbors_list,nan_list,'rows');
  
else
  neighbors_list=[];
end
end
end