function [quiver_handle] = plot_quiver_stacked_scaled(theta, theta2, point_locations, scale, scale2, refPnt)
%PLOT_QUIVER_1 Make scaled quiverplot
%   Detailed explanation goes here
%   point_locations 
%% I. SETUP


% refPnt = [42, 59];
% refPnt = [-35, 35];
% refPnt = ref_point; % jan sigurd's data

% parse reference point
refX = refPnt(1,1); refY = refPnt(1,2);

% we want it to be nx1 (n is # of observations)
theta = theta'; theta2 = theta2';
scale = scale'; scale2 = scale2';

% make sure that point_locations and theta go in the
% same direction (add code for this)

% parse inputs
x = point_locations(:,1);
y = point_locations(:,2);


% calculate the offset angle for each point
offset = rem(atan2d(x-refX, y-refY)+180, 360);
% offset = rem(atan2d(refY-y, refX-x)+180, 360);
% offset = 180-(offset); % look at the other angle (pythag)

% offset = 180-(offset+90); % look at the other angle (pythag)


% compute theta-offset
% theta_shifted = theta-offset;
theta_shifted = theta;
% correct for negative angles if any are found
neg_idx = find(theta_shifted<0);
theta_shifted(neg_idx) = theta_shifted(neg_idx)+360;

% scale the scaling factor
fac = 5;
scale = scale.*fac; scale2 = scale2.*fac;

% Data is organized as (x, y, theta in degrees)
data = [x, y, theta_shifted, scale];
data2 = [x, y, theta2, scale2];

% find [unscaled] vector components (u,v)
% Based on equations: x = x0 + r*cos(theta), y = y0 + r*sin(theta)
u = cos(data(:,3) * pi/180); u2 = cos(data2(:,3) * pi/180);
v = sin(data(:,3) * pi/180); v2 = sin(data2(:,3) * pi/180);

% find the scaling factor (sf)
sf = abs(scale./(sqrt((u.^2)+(v.^2))));
sf2 = abs(scale2./(sqrt((u2.^2)+(v2.^2))));

% multiply components by scaling factor
uprime = u.*sf; uprime2 = u2.*sf2;
vprime = v.*sf; vprime2 = v2.*sf2;

%% II. PLOT
% format figure
quiver_handle = figure;
pbaspect([1 1 1])
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

% plot theta2
h11 = quiver(data2(:,1), data2(:,2), uprime2, vprime2, 0);
set(h11, 'Color', 'k', 'AutoScale', 'off')

hold off;
end

%% SCRATCH CODE

% check if scaling is messing anything up
% h2 = quiver(data(:,1), data(:,2), u, v, 0); 
% set(h2, 'Color', 'b', 'AutoScale', 'on')

