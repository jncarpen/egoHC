function plot_modelfit(out, P, ST, Z)
%   INPUT -
%   Z:          angular variable (HD or MD) ranging from -180 to 180 deg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract behavioral variables
t = P(:,1);
tpf = mode(diff(t)); % time per frame (s)
X = P(:,2);
Y = P(:,3);
Z = rad2deg(Z);

% average cosine (with model parameters)
g = out.model.fitParams.g;
thetaP = out.model.fitParams.thetaP;
xref = out.model.fitParams.xref;
yref = out.model.fitParams.yref;
ang_bins = linspace(-180,180,36);
tc_model = 1 + g*cos(mod(ang_bins-thetaP,360)-180);

% calculate egocentric bearing values (data)
ego_data = mod(atan2d(yref-Y, xref-X) - Z, 360)-180;

% egocentric bearing at time of spikes
idx = knnsearch(P(:,1), ST);
spkX = X(idx); spkY = Y(idx);
spk_ego = ego_data(idx);

% generate tuning curve
edges = linspace(-180, 180, 36);
[spkMap, mapAxis] = histcounts(spk_ego, edges);
[angMap] = histcounts(ego_data, edges);
for i = 1:length(mapAxis)
    if i+1 <= length(mapAxis)
        binCtrs_egoAng(i) = ((mapAxis(i+1)-mapAxis(i))/2)+mapAxis(i);
    end
end
tc_data = (spkMap./(angMap*tpf))';
% tc_data = imgaussfilt(tc_data, 2, 'Padding', 'circular');

figure; hold on;
plot(mapAxis,tc_model');
plot(mapAxis,tc_data, '.');
































% get real modulation curve
for h = 1:10
    real_tc(h) = mean(reshape(out.data(:,:,h), 100, 1));
    real_std(h) = std(reshape(out.data(:,:,h), 100, 1));
end

% make tuning curve (fit by model)
for h = 1:10
    fitted_tc(h) = mean(reshape(out.pred(:,:,h), 100, 1));
end

% calculate the upper & lower bounds of y (+/- 1 STD)
yu = real_tc + 1*real_std;
yl = real_tc - 1*real_std;

% figure

set(gca, 'FontName', 'Calibri Light', 'FontSize', 14);
hold on;
% standard_dev = fill([deg2rad(model.bins-180) fliplr(deg2rad(model.bins-180))], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none');
% alpha(standard_dev, 0.5)

plt_real = plot(deg2rad(out.bins-180), real_tc,'k','LineWidth', .90);
plt_fit = plot(deg2rad(out.bins-180), fitted_tc, '--r', 'LineWidth', .90);

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
th_rad = deg2rad(out.bestParams.thetaP-180);
if th_rad > 0
    plot(t,vertshift + out.bestParams.g*cos(1+t-th_rad), '--r', 'LineWidth', 1.1)
else th_rad < 0;
    plot(t,vertshift + out.bestParams.g*cos(1+t+th_rad), '--r', 'LineWidth', 1.1)
end
hold on;
plt_real = plot(deg2rad(out.bins-180), real_tc,'k','LineWidth', .90);


end
