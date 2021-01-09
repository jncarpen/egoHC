function [err] = cosErr(p,X,Y,H,R)
% Fit a cosine function to a head-direction tuning curve.
%   p:          parameters to fit
%   o:          other parameters
%   R:          data (hd tuning curve)
% @todo need to impose the constraint that pred > 0 **
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model (predicted value)
pred = 1 + (p.g).*cosd((atan2d(p.yref-Y, p.xref-X)+180 - H) - p.thetaP);
pred(pred == 0) = NaN;

% calculate the sums of squared error between 
% the model's prediction and the data
err = nansum((pred(:)-R(:)).^2);


end

%% SCRATCH CODE
% % this is the model proposed in the 'methods' section of Jercog et al.
% pred = 1 + g*(cosd(((y_ref - y)/(x_ref - x)) - theta_P) - (pdx(1)*(cosd(ang(1)-theta_P)) + pdx(2)*(cosd(ang(2)-theta_P)) + ...
%    pdx(3)*(cosd(ang(3)-theta_P)) + pdx(4)*(cosd(ang(4)-theta_P)) + pdx(5)*(cosd(ang(5)-theta_P)) + ...
%    pdx(6)*(cosd(ang(6)-theta_P)) + pdx(7)*(cosd(ang(7)-theta_P)) + pdx(8)*(cosd(ang(8)-theta_P)) + ...
%    pdx(9)*(cosd(ang(9)-theta_P))));


