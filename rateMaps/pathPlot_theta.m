function pathPlot_theta(pos_, SpikeTimes_, rawEEG)
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
if length(pos_) > 3
    pos_ = pos_(:,1:3);
end

t = pos_(:,1);
x = pos_(:,2); % grab xpos
y = pos_(:,3); % grab ypos

% grab theta phase
[~, spkThetaPhz, ~] = getThetaPhz(rawEEG, SpikeTimes_);

% get into degrees 
spkThetaPhz_deg = rad2deg(spkThetaPhz + pi);

spkPos = [];
spkPhz = [];
col1 = [];
col2 = [];
col3 = [];

[spkPos] = data.getSpikePositions(SpikeTimes_', pos_); % spkPos: [t x y]

idx = knnsearch(t, SpikeTimes_);
spkPhz(:,1) = t(idx);
spkPhz(:,2) = spkThetaPhz_deg;

% Adjust vectors to match if needed
if length(spkPos) > length(spkPhz)
    idx = knnsearch(spkPos(:,1), spkPhz(:,1));
    col1 = spkPos(:,1);
    col2 = spkPos(:,2);
    col3 = spkPos(:,2);
    spkPos = [col1(idx), col2(idx), col3(idx)];
elseif length(spkPos) < length(spkPhz)
    idx = knnsearch(spkPhz(:,1), spkPos(:,1));
    col1 = spkPhz(:,1);
    col2 = spkPhz(:,2);
    spkPhz = [col1(idx), col2(idx)];
end

if size(spkPos,2) == 3 && size(spkPhz,2) == 2
else
    disp("Error: spkPos and spkAng have the wrong dimensions.")
end


% make path plot
% figure
plot(x, y, 'Color', [.7 .7 .7])
hold on
scatter(spkPos(:,2), spkPos(:,3), [50], spkPhz(:,2), '.')
% newmap = brighten(hsv,-.8);
colormap(hsv);
colorbar;
caxis([0 360])
pbaspect([1 1 1])
xlim([nanmin(pos_(:,2)), nanmax(pos_(:,2))])
ylim([nanmin(pos_(:,3)), nanmax(pos_(:,3))])
title("Theta phase Path Plot")
box off
end

