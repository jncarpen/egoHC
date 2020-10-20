function [spkSpd] = getSpikeSpeed(pos_,SpikeTimes_)
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

% get speed (in cm/s)
numLeds = 2; 
% parse and medfilt pos vector to eliminate major outliers
t = pos_(:,1); % in seconds
x = medfilt1(pos_(:,2)); 
x2 = medfilt1(pos_(:,4));
y = medfilt1(pos_(:,3)); 
y2 = medfilt1(pos_(:,5));
v = zeros(size(pos_,1), numLeds); % velocity

for i = 2:size(pos_,1)-1
v(i, 1) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
v(i, 2) = sqrt((x2(i+1) - x2(i-1))^2 + (y2(i+1) - y2(i-1))^2) / (t(i+1) - t(i-1));
end

% pad the vector
v(1,1) = v(2,1);
v(1,2) = v(2,2);
v(end, 2) = v(end-1, 2);
v(end, 1) = v(end-1, 1);        
        
% grab speed for first LED
spdLED1 = v(:,1);
minSpd = nanmin(spdLED1);
maxSpd = nanmax(spdLED1);
t = pos_(:,1);
sampleRate = mode(diff(t));

% find closest timestamp & store values
idx = knnsearch(t, SpikeTimes_);
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
    
    % smooth tuning curve values
    tcVals = imgaussfilt(tcVals, 2, 'Padding', 'replicate');


    % plot
    plot(binCtrs, tcVals, 'Color', 'k', 'LineWidth', 1.1)
    title("Speed TC")
    xlabel("speed (cm/s)")
    ylabel("fr (Hz)")
    xlim([0 nanmax(spdLED1)])
    box off
    
    return

end

