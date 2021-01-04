function [head_direction] = get_hd(position)
%GET_HD 
%   Compute head direction (in degrees)
%   INPUTS
%   'position'          position vector with 2 LEDS; looks like this:
%                       [t x y x2 y2]
%   OUTPUTS
%   'head_direction'    head direction (same length as the 'position' input
%                       vector), calculated in degrees.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse position input
posx = position(:,2); 
posy = position(:,3);
posx2 = position(:,4);
posy2 = position(:,5);

% compute head direction in degrees
head_direction = rem(atan2d(posy2-posy,posx2-posx)+180, 360);

end

