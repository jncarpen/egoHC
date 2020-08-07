function [spkSpd] = getSpikeSpeed(pos, speed,SpikeTimes)
%GETSPIKESPEED Find speed of the animal each time the cell fired.
%   INPUT
%   speed:              [spdLED1 spdLED2]
%   SpikeTimes:         1xS vector, where S is number of times the cell fired.
%   
%   OUTPUT
%   spkSpd:             1xS vector (should be same length as SpikeTimes)
%                       with the speed value (in cm/s) of the animal at each timestamp that the
%                       cell fired.
%
%   Jordan Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab speed for first LED
spdLED1 = speed(:,1);
minSpd = nanmin(spdLED1);
maxSpd = nanmax(spdLED1);
t = pos(:,1);
sampleRate = mode(diff(t));

% find closest timestamp & store values
idx = knnsearch(t, SpikeTimes);
spkSpd = spdLED1(idx);

%% calculate speed tuning curve

    nBins = 10;
    edges = linspace(0,300,nBins); % get rid of outliers?
    [spkSpdmap, mapAxis] = histcounts(spkSpd,edges);
    [allSpdmap] = histcounts(spdLED1,edges);
    
    % compute bin centers
    for i = 1:length(mapAxis)
        if i+1 <= length(mapAxis)
            binCtrs(i) = ((mapAxis(i+1)-mapAxis(i))/2)+mapAxis(i);
        end
    end
    
    % calculate tuning curve values
    tcVals = spkSpdmap./(allSpdmap*sampleRate + eps); 
    
    % plot
    plot(binCtrs, tcVals, 'Color', 'k', 'LineWidth', 1.5)
    title("Speed TC")
    xlabel("speed (units/s)")
    ylabel("fr (Hz)")
    xlim([0 300])
    % xticks([-pi -pi/2 0 pi/2 pi])
    % xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
    box off
    
    return

end

