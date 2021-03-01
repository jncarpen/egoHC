function plot_vonMises(theta)
%PLOT_VONMISES Make a sample plot of von-Mises distribution

% default parameters
alpha = 0:0.001:2*pi; % from 0 to 2pi
theta = deg2rad(theta);
kappa = 2;

% calculate the pdf
[p, ~] = circ_vmpdf(alpha, theta, kappa);

% plot
figure
hold on;
set(gcf,'color','w');
plot(rad2deg(alpha), p, 'Color','k', 'LineWidth', 1.5);
xline(rad2deg(theta), '--r');
pbaspect([1 1 1])
xlim([0 360])
xticks([0 90 180 270])
xlabel("angle (degrees)", 'FontName', 'Calibri light', 'FontSize', 12, 'FontWeight', 'normal')
ylabel("pdx", 'FontName', 'Calibri light', 'FontSize', 12, 'FontWeight', 'normal')
title("von-mises distribution", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
box off
hold off;

end

