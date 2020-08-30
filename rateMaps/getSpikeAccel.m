function [spkAcc] = getSpikeAccel(pos, accel, SpikeTimes)
%GETSPIKESPEED Find speed of the animal each time the cell fired.
%   INPUT
%   accel:              [accLED1 accLED2]
%   SpikeTimes:         1xS vector, where S is number of times the cell fired.
%   
%   OUTPUT
%   spkAccel:           1xS vector (should be same length as SpikeTimes)
%                       with the acceleration value (in cm/s) of the animal at 
%                       each timestamp that the cell fired.
%                       
%
%   Jordan Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"))

% grab speed for first LED
accLED1 = accel(:,1);
minAcc = nanmin(accLED1);
maxAcc = nanmax(accLED1);
t = pos(:,1);
sampleRate = mode(diff(t));

% find closest timestamp & store values
idx = knnsearch(t, SpikeTimes);
spkAcc = accLED1(idx);

%% calculate acceleration tuning curve

    nBins = 10;
    edges = linspace(0,300,nBins); % get rid of outliers?
    [spkAccmap, mapAxis] = histcounts(spkAcc,edges);
    [allAccmap] = histcounts(accLED1,edges);
    
    % compute bin centers
    for i = 1:length(mapAxis)
        if i+1 <= length(mapAxis)
            binCtrs(i) = ((mapAxis(i+1)-mapAxis(i))/2)+mapAxis(i);
        end
    end
    
    % calculate tuning curve values
    tcVals = spkAccmap./(allAccmap*sampleRate + eps); 
    % smooth_tcVals = general.smooth(tcVals, [2 2]);
    
    % plot
    plot(binCtrs, tcVals, 'Color', 'k', 'LineWidth', 1.5)
    % plot(binCtrs, tcVals, 'Color', 'k', 'LineWidth', 1.5)
    title("Acc TC")
    xlabel("acc (units/s)")
    ylabel("fr (Hz)")
    % xlim([0 300])
    % xticks([-pi -pi/2 0 pi/2 pi])
    % xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
    box off
    
    return

end

