function plot_cosineFit(out, P, ST)
%PLOT_COSINEFIT

% best fit parameters
g = out.model.fitParams.g;
thetaP = out.model.fitParams.thetaP;
xref = out.model.fitParams.xref;
yref = out.model.fitParams.yref;

best = out.bestParams;
X = out.loc3d(:,:,1); 
Y = out.loc3d(:,:,2);

% make the head direction array
nBins = 10;
H = zeros(nBins, nBins, nBins);
for h = 1:length(out.bins)
    H(:,:,h) = ones(10,10).*out.bins(h);
end

% calculate the cosine curve
cosine = (best.g).*cosd(1+(atan2d(best.yref-Y, best.xref-X)+180 - H) - best.thetaP);

% tuning curve of actual data (HD)
for h = 1:10
    real_tc(h) = nanmean(reshape(out.data(:,:,h), 100, 1));
    real_std(h) = std(reshape(out.data(:,:,h), 100, 1));
end

% calculate the upper & lower bounds of y (+/- 2 stds)
yu = real_tc + 1*real_std;
yl = real_tc - 1*real_std;

hold on;
standard_dev = fill([deg2rad(out.bins-180) fliplr(deg2rad(out.bins-180))], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none');
alpha(standard_dev, 0.5)
set(gca, 'FontName', 'Calibri Light', 'FontSize', 14);
plt_real = plot(deg2rad(out.bins-180), real_tc,'ko','MarkerFaceColor','b');
plt_fit = plot(deg2rad(out.bins-180),fitted_tc, '--r', 'LineWidth', .90);
xlabel("egocentric angle (rad)", "Fontsize", 14)
ylabel("activity (normalized)", "Fontsize", 14)
title("Egocentric Response", 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
% l.FontSize = 12; l.FontName = 'Calibri light'; l.Location = 'northeastoutside';
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
pbaspect([1 1 1])
box off

%% FUNCTIONS
function [head_direction] = get_hd(position)
    posx = position(:,2); posy = position(:,3);
    posx2 = position(:,4); posy2 = position(:,5);
    % compute head direction in degrees
    head_direction = rem(atan2d(posy2-posy,posx2-posx)+180, 360);
end
end

