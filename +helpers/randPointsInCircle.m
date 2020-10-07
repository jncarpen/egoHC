function [x,y] = randPointsInCircle(x1, y1, rc)
%RANDPOINTSINCIRCLE Summary of this function goes here
%   Generate n random points within the radius of a circle
%   INPUTS
%   'rc'        radius constant
%   'x1'         center x
%   'y1'         center y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=2*pi*rand;
r=sqrt(rand);
x=(rc*r)*cos(a)+x1;
y=(rc*r)*sin(a)+y1;

end
