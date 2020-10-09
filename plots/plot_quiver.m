function plot_quiver(theta, point_locations)
%PLOT_QUIVER Summary of this function goes here
%   Detailed explanation goes here

x = point_locations(:,1);
y = point_locations(:,2);
theta = theta';

figure
pbaspect([1 1 1])
box off

% Data is organized as (x, y, theta in degrees)
data = [x, y, theta];

% Identify data from 1st quadrant and store in q1, identify data from 2nd quadrant
% goes to q2, and so on. A value of 1 indicates it is in the quadrant, 0 otherwise
q1 = (data(:,3) <= 90);
q2 = (data(:,3) > 90) .* (data(:,3) <= 180);
q3 = (data(:,3) > 180) .* (data(:,3) <= 270);
q4 = (data(:,3) > 271) .* (data(:,3) <= 360);
hold on;

% Use QUIVER to specify the start point (tail of the arrow) and direction based on angle
% q1, q2, q3, and q4 are used to generate four different QUIVER handles (h1, h2, h3, and h4)
% This is necessary for varying colors based on direction
% Based on equations: x = x0 + r*cos(theta), y = y0 + r*sin(theta)
% In the usage below, x0 = data(:,1), y0 = data(:,2), theta = data(:,3) * pi / 180
% Can also specify a scale factor as the last argument to quiver (not specified below)
h1 = quiver(data(q1 == 1,1), data(q1 == 1,2), cos(data(q1 == 1,3) * pi/180), sin(data(q1 == 1,3) * pi/180));
h2 = quiver(data(q2 == 1,1), data(q2 == 1,2), cos(data(q2 == 1,3) * pi/180), sin(data(q2 == 1,3) * pi/180)); % sin is negative in 2nd quadrant
h3 = quiver(data(q3 == 1,1), data(q3 == 1,2), cos(data(q3 == 1,3) * pi/180), sin(data(q3 == 1,3) * pi/180));
h4 = quiver(data(q4 == 1,1), data(q4 == 1,2), cos(data(q4 == 1,3) * pi/180), sin(data(q4 == 1,3) * pi/180)); % cos is negative in 4th quadrant

% Set colors to red for 1st quadrant, blue for 2nd, green for 3rd, cyan for 4th
% Also, turn scaling off. get(h1) will return additional property-value pairs
% set(h1, 'Color', [.8 0 .5], 'AutoScale', 'off')
% set(h2, 'Color', [.9 .5 .3], 'AutoScale', 'off')
% set(h3, 'Color', [.1 0 .8], 'AutoScale', 'off')
% set(h4, 'Color', [0 .75 .5], 'AutoScale', 'off')


%
set(h1, 'Color', 'k', 'AutoScale', 'on')
set(h2, 'Color', 'k', 'AutoScale', 'on')
set(h3, 'Color', 'k', 'AutoScale', 'on')
set(h4, 'Color', 'k', 'AutoScale', 'on')

% Done plotting
hold off;

end

