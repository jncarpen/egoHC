function [simulatedSpikes,simulatedSpk2PosInd] = spatialPoissonSimulation(P,spikes)
%spatialPoissonSimulation
%   Generate a simulated spike train based on the expected firing of the
%   cell at each location (rate map). This is the method used by Wang et
%   al. (2018) to test how much egocentric tuning we expect to find only
%   based on the animal's trajectory and the cell's spatial firing pattern.
%   Inputs: <pos> Nx3 matrix of position data (timestamps, x, y)
%           <spikes> Mx3 matrix of spike data (timestamps, x, y), where M<N
%   Note: to get HD at time of simulated spike timestamps, just write
%   HD(simulatedSpk2PosInd) where <HD> is an Nx1 vector (same length as
%   position data) with animal's head direction


%calculate sample time
sampleTime=mode(diff(P(:,1))); 

%calculate rate map in spikes/second
map=analyses.map(P,spikes(:,1));

%multiply by sample time to get spikes/frame (or spikes/position sample)
lambdaMatrix=map.z.*sampleTime; 

%get x- and y bins 
xNumBins=size(map.z,2);
yNumBins=size(map.z,1); 

%prepare
simulatedSpikes=[]; 
simulatedSpk2PosInd=[]; 

%get number of samples
numSamples=size(P,1); 

%iterate over position samples
for z=1:numSamples
    
    %get current position sample
    currentSample=P(z,:); 
    
    %find x bin of position sample
    for j=1:xNumBins
        
        if currentSample(2)>=map.x(j) && currentSample(2)<map.x(j+1)
            xBinInd=j;
            break
        else
            xBinInd = 0;
        end 
        
    end
    
    %find y bin of position sample  
    for i=1:yNumBins
        
        if currentSample(3)>=map.y(i) && currentSample(3)<map.y(i+1)
            yBinInd=i;
            break
        else
            yBinInd = 0;
        end 
        
    end 
    
    if yBinInd > 0 && xBinInd > 0
        %get lambda (firing rate in frames/second) corresponding to current bin 
        currentLambda=lambdaMatrix(yBinInd,xBinInd); 

        %draw random number of spikes (n) from Poisson distribution
        n=poissrnd(currentLambda); 

        %store number of spikes and spike-to-position indices
        if n>0
            simulatedSpikes=[simulatedSpikes;repmat(currentSample,n,1)];
            simulatedSpk2PosInd=[simulatedSpk2PosInd;repmat(z,n,1)]; 
        end 
    end
    
end 

end


%% make a ratemap

% x/y pos
xpos = P(:,2);
ypos = P(:,3);

% bin the arena and find bin centers
nBins = 100;
xEdges = linspace(0,150,nBins+1);
yEdges = linspace(0,150,nBins+1);
xCenter = (diff(xEdges)./2)+xEdges(1:end-1);
yCenter = (diff(yEdges)./2)+yEdges(1:end-1);
[X,Y] = meshgrid(xCenter,yCenter);

% get variance
xvar = nanvar(xpos); 
yvar = nanvar(ypos); 

% get standard deviation (sigma)
xsd = sqrt(xvar); 
ysd = sqrt(yvar);

% kendall or pearson correlation (excludes NaNs)
warning('off','all')
rho = corr(xpos,ypos,'Type','Pearson','Rows','complete');
 if (abs(rho) >= 1.0)
        disp("error: rho must lie between -1 and 1");
    return
 end
 warning('on','all')
 
% % calculation of the covariance matrix
% covxy = rho*xsd*ysd;
% scaleCov = .75; % scale cov by some factor
% C = [xvar covxy; covxy yvar].*scaleCov; % the covariance matrix
% A = inv(C); % the inverse covariance matrix

w = 1000;
C = [w 0; 0 w];
A = inv(C);

xmu = 10; ymu = 10;


% Compute value of Gaussian pdf at each point in the grid
clear zMapped
% z = (1/sqrt(2*pi*det(C))) * exp(-1/2* (A(1,1)*(X-xmu).^2 + 2*A(1,2)*(X-xmu).*(Y-ymu) + A(2,2)*(Y-ymu).^2));
z = 10.* exp(-1/2* (A(1,1)*(X-xmu).^2 + 2*A(1,2)*(X-xmu).*(Y-ymu) + A(2,2)*(Y-ymu).^2));

figure
imagesc(z); colorbar;































%%
% A: amplitude
% x0, y0: center 
% sigma_x, sigma_y: x and y spreads of the blob
y0 = 10; x0=10;
sigmaX=5; sigmaY=10;
A = 10; % firing rate (Hz)
fxy = A.*exp(-(((X-x0).^2)./(2*sigmaX.^2) + ((Y-y0).^2.)/(2*sigmaY.^2)));

figure
surf(fxy)
















