function plotThetaPhzTC(tc_TP)
%PLOTTHETAPHZTC Summary of this function goes here
%   INPUTS
%   tc_TP:              theta phase tuning curve, 2xB vector, where B is the
%                       number of bins. This variable can be generated via the
%                       getThetaPhz.m function.
%
%   OUTPUTS
%   tuningCurve:        Theta phase tuning curve will be generated with this function.
%
%   Jordan Carpenter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse tc_TP vector
binCtrs = tc_TP(:,1);
tcVals = tc_TP(:,2);

% smooth tcVals
% tcVals = smoothdata(tcVals, 'gaussian', 10);
tcVals = imgaussfilt(tcVals, 2, 'Padding', 'circular');

% plot
plot(binCtrs, tcVals, 'Color', 'k', 'LineWidth', 1.10)
title("Theta Phase TC")
xlabel("phase (rad)")
ylabel("fr (Hz)")
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
box off
end

