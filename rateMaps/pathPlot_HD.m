function pathPlot_hd(pos, SpikeTimes, hd)
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
col1 = [];
col2 = [];
col3 = [];

[spkPos] = data.getSpikePositions(SpikeTimes', pos); % spkPos: [t x y]

idx = knnsearch(t, SpikeTimes);
spkAng(:,1) = t(idx);
spkAng(:,2) = hd(idx);

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
plot(x, y, 'Color', [.85 .85 .85]);
hold on;
scatter(spkPos(:,2), spkPos(:,3), [45], spkAng(:,2), '.')
% newmap = brighten(hsv,-.8);
colormap(hsv);
c = colorbar; c.FontSize = 15; c.Ticks = [90 180 270]; c.Box = 'off';
caxis([0 360]);
c.FontName = 'Helvetica UI'; c.FontWeight = 'normal';
pbaspect([1 1 1])
xlim([nanmin(pos(:,2)), nanmax(pos(:,2))])
ylim([nanmin(pos(:,3)), nanmax(pos(:,3))])
% title("spike plot (colored by HD)", 'FontName', 'Calibri light', 'FontSize', 18, 'FontWeight', 'normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
box off
% hold off;
end

