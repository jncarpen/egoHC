function [quiver_handle] = plot_quiver_stacked(theta, theta2, point_locations)
%PLOT_QUIVER_1 Summary of this function goes here
%   Detailed explanation goes here

% parse inputs
x = point_locations(:,1);
y = point_locations(:,2);
theta = theta'; theta2 = theta2';

% format figure
quiver_handle = figure;
% title("blue:theta, black:theta2")
title("blue:peak(ego), black:peak(hd)")
pbaspect([1 1 1])
box off

% Data is organized as (x, y, theta in degrees)
data = [x, y, theta];
data2 = [x, y, theta2];

% Identify data from 1st quadrant and store in q1, identify data from 2nd quadrant
% goes to q2, and so on. A value of 1 indicates it is in the quadrant, 0 otherwise
q1 = (data(:,3) <= 90);
q2 = (data(:,3) > 90) .* (data(:,3) <= 180);
q3 = (data(:,3) > 180) .* (data(:,3) <= 270);
q4 = (data(:,3) > 271) .* (data(:,3) <= 360);

q11 = (data2(:,3) <= 90);
q22 = (data2(:,3) > 90) .* (data2(:,3) <= 180);
q33 = (data2(:,3) > 180) .* (data2(:,3) <= 270);
q44 = (data2(:,3) > 271) .* (data2(:,3) <= 360);

hold on;

% Use QUIVER to specify the start point (tail of the arrow) and direction based on angle
% q1, q2, q3, and q4 are used to generate four different QUIVER handles (h1, h2, h3, and h4)
% This is necessary for varying colors based on direction
% Based on equations: x = x0 + r*cos(theta), y = y0 + r*sin(theta)
% In the usage below, x0 = data(:,1), y0 = data(:,2), theta = data(:,3) * pi / 180
% Can also specify a scale factor as the last argument to quiver (not specified below)
h1 = quiver(data(q1 == 1,1), data(q1 == 1,2), cos(data(q1 == 1,3) * pi/180), sin(data(q1 == 1,3) * pi/180), 'color',[1 0 1]);
h2 = quiver(data(q2 == 1,1), data(q2 == 1,2), cos(data(q2 == 1,3) * pi/180), sin(data(q2 == 1,3) * pi/180), 'color',[1 0 1]); % sin is negative in 2nd quadrant
h3 = quiver(data(q3 == 1,1), data(q3 == 1,2), cos(data(q3 == 1,3) * pi/180), sin(data(q3 == 1,3) * pi/180), 'color',[1 0 1]);
h4 = quiver(data(q4 == 1,1), data(q4 == 1,2), cos(data(q4 == 1,3) * pi/180), sin(data(q4 == 1,3) * pi/180), 'color',[1 .25 1]); % cos is negative in 4th quadrant

h11 = quiver(data2(q11 == 1,1), data2(q11 == 1,2), cos(data2(q11 == 1,3) * pi/180), sin(data2(q11 == 1,3) * pi/180), 'color',[0 0 0]);
h22 = quiver(data2(q22 == 1,1), data2(q22 == 1,2), cos(data2(q22 == 1,3) * pi/180), sin(data2(q22 == 1,3) * pi/180), 'color',[0 0 0]); % sin is negative in 2nd quadrant
h33 = quiver(data2(q33 == 1,1), data2(q33 == 1,2), cos(data2(q33 == 1,3) * pi/180), sin(data2(q33 == 1,3) * pi/180), 'color',[0 0 0]);
h44 = quiver(data2(q44 == 1,1), data2(q44 == 1,2), cos(data2(q44 == 1,3) * pi/180), sin(data2(q44 == 1,3) * pi/180), 'color',[0 0 0]); % cos is negative in 4th quadrant

% set properties for theta
set(h1, 'Color', 'k', 'AutoScale', 'on')
set(h2, 'Color', 'k', 'AutoScale', 'on')
set(h3, 'Color', 'k', 'AutoScale', 'on')
set(h4, 'Color', 'k', 'AutoScale', 'on')

% set properties for theta2
set(h1, 'Color', 'b', 'AutoScale', 'on')
set(h2, 'Color', 'b', 'AutoScale', 'on')
set(h3, 'Color', 'b', 'AutoScale', 'on')
set(h4, 'Color', 'b', 'AutoScale', 'on')

% Done plotting
% hold off;

end

