function plot_modelfit(model)

% get real modulation curve
for h = 1:10
    real_tc(h) = mean(reshape(model.data(:,:,h), 100, 1));
    real_std(h) = std(reshape(model.data(:,:,h), 100, 1));
end

% make tuning curve (fit by model)
for h = 1:10
    fitted_tc(h) = mean(reshape(model.pred(:,:,h), 100, 1));
end

% calculate the upper & lower bounds of y (+/- 1 STD)
yu = real_tc + 1*real_std;
yl = real_tc - 1*real_std;

% figure

set(gca, 'FontName', 'Calibri Light', 'FontSize', 14);
hold on;
% standard_dev = fill([deg2rad(model.bins-180) fliplr(deg2rad(model.bins-180))], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none');
% alpha(standard_dev, 0.5)

plt_real = plot(deg2rad(model.bins-180), real_tc,'k','LineWidth', .90);
plt_fit = plot(deg2rad(model.bins-180), fitted_tc, '--r', 'LineWidth', .90);

xlabel("egocentric angle (rad)", "Fontsize", 14)
ylabel("activity (normalized)", "Fontsize", 14)
title("Model Fit", 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
% l.FontSize = 12; l.FontName = 'Calibri light'; l.Location = 'northeastoutside';
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
pbaspect([1 1 1])
box off

figure
t = -pi:.01:pi;
vertshift = .4;
th_rad = deg2rad(model.bestParams.thetaP-180);
if th_rad > 0
    plot(t,vertshift + model.bestParams.g*cos(1+t-th_rad), '--r', 'LineWidth', 1.1)
else th_rad < 0;
    plot(t,vertshift + model.bestParams.g*cos(1+t+th_rad), '--r', 'LineWidth', 1.1)
end
hold on;
plt_real = plot(deg2rad(model.bins-180), real_tc,'k','LineWidth', .90);


end
