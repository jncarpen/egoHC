function [C] = circ_add(A,B)
%CIRC_SUBTRACT Subtract circular varibles (deg)
%   A+B (vectors)

C = mod(A+B, 360);
% C = C(C<0);
% C = C(C<0) + (round(abs(C((C<0))/360)).*360); % neg values

end


