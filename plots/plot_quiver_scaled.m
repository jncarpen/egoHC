function [quiver_handle] = plot_quiver_scaled(theta, point_locations, scale, refPnt, offset)
%PLOT_QUIVER_1 Make scaled quiverplot
%   Detailed explanation goes here
%   point_locations 
%% I. SETUP

% parse reference point
refX = refPnt(1,1); refY = refPnt(1,2);

% we want it to be nx1 (n is # of observations)
theta = theta'; 
scale = scale'; 

% make sure that point_locations and theta go in the
% same direction (add code for this)

% parse inputs
x = point_locations(:,1);
y = point_locations(:,2);

switch offset
    case "True"
        % calculate the offset angle for each point
        offset = rem(atan2d(x-refX, y-refY)+180, 360);
        % offset = rem(atan2d(refY-y, refX-x)+180, 360);
        offset = 180-(offset); % look at the other angle (pythag)
        % offset = 180-(offset+90); % look at the other angle (pythag)
        theta_shifted = theta-offset;
    case "False"
        theta_shifted = theta;

% correct for negative values
neg_idx = find(theta_shifted<0);
theta_shifted(neg_idx) = theta_shifted(neg_idx)+360;

% scale the scaling factor
fac = 5;
scale = scale.*fac;

% Data is organized as (x, y, theta in degrees)
data = [x, y, theta_shifted, scale];

% find [unscaled] vector components (u,v)
% Based on equations: x = x0 + r*cos(theta), y = y0 + r*sin(theta)
u = cos(data(:,3) * pi/180); 
v = sin(data(:,3) * pi/180); 

% find the scaling factor (sf)
sf = abs(scale./(sqrt((u.^2)+(v.^2))));

% multiply components by scaling factor
uprime = u.*sf; 
vprime = v.*sf; 

%% II. PLOT
% format figure
quiver_handle = figure;
% pbaspect([1 1 1])
box off
hold on;

% plot theta 1 
h1 = quiver(data(:,1), data(:,2), uprime, vprime, 0); % 0 turns autoscaling off
set(h1, 'Color', 'k', 'AutoScale', 'off')

% plot reference location
refPnt_plot = plot(refPnt(1,1), refPnt(1,2), 'o', 'MarkerSize', 10);
set(refPnt_plot, 'markerfacecolor', 'blue');
boxSize = 150;
xlim([0 boxSize])
ylim([0 boxSize])

hold off;
end
