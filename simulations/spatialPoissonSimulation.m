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

