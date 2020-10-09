function pathPlot(pos, SpikeTimes)
%PATHPLOT: make simple pathplot

%   INPUT:
%   pos:        [t x y] or [t x y x2 y2]
%   spikeTs:    Timestamps (in seconds) for a single neuron

%   OUTPUT:
%   pathPlot:   Function will return pathplot with spikes overlaid

%   Jordan Carpenter, 2020.

%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%

% add BNT to path
addpath(genpath('C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001'));

% remove [x2 y2] if present
if length(pos) > 3
    pos = pos(:,1:3);
end

x = pos(:,2); % grab xpos
y = pos(:,3); % grab ypos

[spkPos] = data.getSpikePositions(SpikeTimes', pos); % spkPos: [t x y]

% make path plot
% figure
plot(x, y, 'Color', [.7 .7 .7])
hold on
scatter(spkPos(:,2),spkPos(:,3), 3, [1, 0, 0], 'filled')
title("Path Plot")
box off
% xlabel("X")
% ylabel("Y")

end

