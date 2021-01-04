function [boxCtrX,boxCtrY] = getBoxCenter(pos_cm)
%BOX_CENTER Summary of this function goes here
%   Get location of box center for a single session

xMin = nanmin(pos_cm(:,2));
xMax = nanmax(pos_cm(:,2));
yMin = nanmin(pos_cm(:,3));
yMax = nanmax(pos_cm(:,2));

boxCtrX = (xMin+xMax)/2;
boxCtrY = (yMin+yMax)/2;
end

%% SCRATCH CODE
% plot vertices of box
% plot(xMin,yMin, 'o')
% % hold on
% plot(xMin,yMax, 'o')
% plot(xMax,yMin, 'o')
% plot(xMax,yMax, 'o')
