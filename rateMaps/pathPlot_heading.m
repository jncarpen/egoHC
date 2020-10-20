function pathPlot_heading(pos, SpikeTimes)
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

% parse position
t = pos(:,1);
x = pos(:,2); % grab xpos
y = pos(:,3); % grab ypos

% grab heading direction values
[md_1, md_2, md_3] = get_moving_direction(pos);
heading_dir = md_1';

spkPos = [];
spkAng = [];
col1 = [];
col2 = [];
col3 = [];

[spkPos] = data.getSpikePositions(SpikeTimes', pos); % spkPos: [t x y]

idx = knnsearch(t, SpikeTimes);
spkAng(:,1) = t(idx);
spkAng(:,2) = heading_dir(idx);

% Adjust vectors to match if needed
if length(spkPos) > length(spkAng)
    idx = knnsearch(spkPos(:,1), spkAng(:,1));
    col1 = spkPos(:,1);
    col2 = spkPos(:,2);
    col3 = spkPos(:,2);
    spkPos = [col1(idx), col2(idx), col3(idx)];
elseif length(spkPos) < length(spkAng)
    idx = knnsearch(spkAng(:,1), spkPos(:,1));
    col1 = spkAng(:,1);
    col2 = spkAng(:,2);
    spkAng = [col1(idx), col2(idx)];
end

if size(spkPos,2) == 3 && size(spkAng,2) == 2
else
    disp("Error: spkPos and spkAng have the wrong dimensions.")
end


% make path plot
% figure
plot(x, y, 'Color', [.7 .7 .7]);
hold on
scatter(spkPos(:,2), spkPos(:,3), [50], spkAng(:,2), '.')
% newmap = brighten(hsv,-.8);
colormap(hsv);
colorbar;
caxis([0 360])
pbaspect([1 1 1])
xlim([nanmin(pos(:,2)), nanmax(pos(:,2))])
ylim([nanmin(pos(:,3)), nanmax(pos(:,3))])
title("Heading Path Plot")
set(gca,'xtick',[])
set(gca,'ytick',[])
box off
end

