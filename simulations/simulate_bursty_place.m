function [sim] = simulate_bursty_place(map, P)
%SIMULATE_PLACE_V2 
%   INPUTS -
%   P:                      position vector (in cm) [t x y]
%   map.z:                  ratemap (Hz)
%   map.whichBin.x:         bin the animal occupies for each frame (x)
%   map.whichBin.y:         bin the animal occupies for each frame (y)
%
%   OUTPUTS -
%   sim.ST:                 simulated spiketimes (s)
%   sim.TR:                 simulated spiketrain
%   J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ratemap
RM = map.z;

% sampling freq (s)
t = P(:,1);
fs = mode(diff(t));

% multiply by sample time to get spikes/frame (or spikes/position sample)
lambdaMatrix = RM'.*fs;

simTrn = zeros(length(t), 1);
simTimes = [];
nSpksBrst = 3; % number of spikes/burst
nFramesBrst = 10; % number of frames/burst
count = 1;
for frame = 1:length(t)
    if count < length(t)
        
        % timestamp now
        tnow = t(count);

        % find x,y bin for current sample
        xnow = map.whichBin.x(count);
        ynow = map.whichBin.y(count);

        % if timestamp has a nan
        if isnan(xnow) || isnan(ynow)
            simTrn(count) = 0;
        else

        % find value of ratemap (fr now)
        lambda = lambdaMatrix(xnow, ynow);

        %draw random number of spikes (n) from Poisson distribution
        n = poissrnd(lambda);

        if count < length(t)-nFramesBrst && n >= 3
            % burst
            simTrn(count:count+nFramesBrst-1) = nSpksBrst;
            simTimes = [simTimes; repmat(tnow,nSpksBrst,1); repmat(t(count+1),nSpksBrst,1);...
                repmat(t(count+2),nSpksBrst,1); repmat(t(count+3),nSpksBrst,1);...
                repmat(t(count+4),nSpksBrst,1)];
            count = count + 4;
        end
        
        % for the case that there are 1 or 2 spikes
        switch n
            case 1
                simTrn(count) = n;
                simTimes = [simTimes; tnow];
                count = count + 1;
            case 2
                simTrn(count) = n;
                simTimes = [simTimes; tnow; tnow];
                count = count + 1;
            otherwise
                count = count +1;
        end
        end
    end 
end


% save variables in struct
sim.TR = simTrn;
sim.ST = simTimes;

end

