function [corrPos] = correct_pos_general(raw_position, boxSize)
%CORRECT_POS_GENERAL Convert position to cm & shift to start at origin
%(0,0).
%   As of now, the box needs to be square and the boxSize input should
%   be a scalar!
%   J. Carpenter- Oct. 12, 2020.

% get length of the box 
xLen = nanmax(raw_position(:,2))-nanmin(raw_position(:,2));
xLen2 = nanmax(raw_position(:,4))-nanmin(raw_position(:,4));
yLen = nanmax(raw_position(:,3))-nanmin(raw_position(:,3));
yLen2 = nanmax(raw_position(:,5))-nanmin(raw_position(:,5));

% calculate a conversion factor
conFac_x = boxSize/xLen; conFac_y = boxSize/yLen; 
conFac_x2 = boxSize/xLen2; conFac_y2 = boxSize/yLen2;

% define time
t = raw_position(:,1); 

% correct position by shifting everything over (to start at ~0)
x = (raw_position(:,2)-nanmin(raw_position(:,2)))*conFac_x; y = (raw_position(:,3)-nanmin(raw_position(:,3)))*conFac_y;
x2 = (raw_position(:,4)-nanmin(raw_position(:,4)))*conFac_x2; y2 = (raw_position(:,5)-nanmin(raw_position(:,5)))*conFac_y2;

% corrected position
corrPos = [t, x, y, x2, y2];

end

