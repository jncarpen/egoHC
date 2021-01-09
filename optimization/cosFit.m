function [pred, err] = cosFit(p,X,Y,H,R)
% Fit a cosine function to a head-direction tuning curve.
%   p:          parameters to fit
%   o:          other parameters
%   R:          data (hd tuning curve)
% @todo need to impose the constraint that pred > 0 **
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model (predicted value)
pred = 1+(p.g).*cosd((atan2d(p.yref-Y, p.xref-X)+180 - H) - p.thetaP);
% pred(pred == 0) = NaN;

% calculate the sums of squared error between 
% the model's prediction and the data
err = nansum((pred(:)-R(:)).^2);


end

% scratch
% pred = abs(1 + (p.g).*cosd(mod(mod(atan2d(p.yref-Y, p.xref-X)+180 - H, 360) - p.thetaP, 360)));
