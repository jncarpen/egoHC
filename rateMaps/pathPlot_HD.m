function pathPlot_HD(pos, SpikeTimes, hd)
%PATHPLOT_HD: Make pathplot with spikes colored by head direction.

%   INPUT:
%   pos:                [t x y] or [t x y x2 y2]
%   SpikeTimes:         Timestamps (in seconds) for a single neuron
%   hd:                 head direction vector (in degrees); should be same
%                       length as pos vector.

%   OUTPUT:
%   pathPlot:           Function will return pathplot with spikes overlaid

%   Jordan Carpenter, 2020.

%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%

% add BNT to path
addpath(genpath('C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001'));

% remove [x2 y2] if present
if length(pos) > 3
    pos = pos(:,1:3);
end

t = pos(:,1);
x = pos(:,2); % grab xpos
y = pos(:,3); % grab ypos

spkPos = [];
spkAng = [];

[spkPos] = data.getSpikePositions(SpikeTimes', pos); % spkPos: [t x y]
idx = knnsearch(t, SpikeTimes);
spkAng(:,1) = t(idx);
spkAng(:,2) = hd(idx);

% Adjust vectors to match if needed
if length(spkPos) > length(spkAng)
    idx = knnsearch(spkPos(:,1), spkAng(:,1));
    spkPos = spkPos(idx);
elseif length(spkPos) < length(spkAng)
    idx = knnsearch(spkAng(:,1), spkPos(:,1));
    spkAng = spkAng(idx);
end


% make path plot
% figure
plot(x, y, 'Color', [.7 .7 .7])
hold on
scatter(spkPos(:,2), spkPos(:,3), [3], spkAng(:,2), '.')
title("HD Path Plot")
colorbar(gca, "hsv", 'eastoutside')
xlabel("X")
ylabel("Y")
return

end

