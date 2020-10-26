function [P_logical, ST_thresh] = speed_threshold(position, ST, threshold)
%SPEED_THRESHOLD Summary of this function goes here
%   return spike times and position values for a given speed threshold
%   'P_logical': gives 

% position
x = position(:,2);
y = position(:,3);

% get speed values for each timestamp
[s, ~] = get_speed(position);

% speed for LED1
s_LED1 = s(:,1); 

% force NaN values for timestamps where speed 
% is below the specified threshold
P_logical = s_LED1 > threshold;
% position(s_idx, 2:5) = NaN;
% P_thresh = position;

t = position(:,1);
minSpd = nanmin(s_LED1);
maxSpd = nanmax(s_LED1);
sampleRate = mode(diff(t));

% find closest timestamp & store values
idx = knnsearch(t, ST);
spkSpd = s_LED1(idx);
% above_thresh = find(spkSpd>threshold);
ST_thresh = ST(spkSpd>threshold);
end

