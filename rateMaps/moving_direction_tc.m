function moving_direction_tc(position, STimes)
%MOVING_DIRECTION_TC Summary of this function goes here
%   Detailed explanation goes here

% compute sampleTime
t = position(:,1);
sampleTime = mode(diff(t));

% get moving direction (3 diff measures)
[md_1, md_2, md_3] = get_moving_direction(position);

% find moving direction of animal when cell spikes
idx = knnsearch(t, STimes);
spk_md1 = md_1(idx);
spk_md2 = md_2(idx);
spk_md3 = md_3(idx);

% calculate tuning curves
smoothFactor = 3; binWidthFactor = 3;
tc_md1 = analyses.turningCurve(spk_md1, md_1, sampleTime, 'smooth', smoothFactor, 'binWidth', binWidthFactor);
tc_md2 = analyses.turningCurve(spk_md2, md_2, sampleTime, 'smooth', smoothFactor, 'binWidth', binWidthFactor);
tc_md3 = analyses.turningCurve(spk_md3, md_3, sampleTime, 'smooth', smoothFactor, 'binWidth', binWidthFactor);

% plot tuning curves
% figure
% set(gcf,'color','w');
hold on;
plot(tc_md1(:,1), tc_md1(:,2), 'Color', 'r', 'LineWidth', 1.10)
plot(tc_md2(:,1), tc_md2(:,2), 'Color', 'k', 'LineStyle',':', 'LineWidth', 1.10)
plot(tc_md3(:,1), tc_md3(:,2), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.10)
xlim([0 360])
% legend('method1', 'method2', 'method3')
title("Movement Direction")
ylabel("fr (Hz)")
xlabel("movement direction")
xticks([0 90 180 270 360])
box off
hold off;

end

