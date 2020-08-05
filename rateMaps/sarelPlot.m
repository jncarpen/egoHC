function [outputArg1,outputArg2] = sarelPlot(pos, hwLoc)
%SARELPLOT Summary of this function goes here
%   INPUTS
%   pos:        [t x y x2 y2]
%   hwLoc:      A single cell from the 1xS cell array called hwLoc, where 
%               S is the number of sessions for a particular animal (hwLoc
%               can be generated from the getWellLoc.m function). This will
%               be a scalar value of the location of the home well for FM
%               trials. Numbers will be in FM coordinates and *not* XY
%               coordinates, for example: 36 or 37.
%   OUTPUTS
%
%
%
%   Jordan Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse position vector
t = pos(:,1);
x = pos(:,2);
y = pos(:,3);



end

