function plot_circ_hist(sDist,nBins)
%PLOT_CIRC_HIST 
%   Plot a circular histogram (using the CircHist function)
%   Frederick Zittrell (2020). CircHist - circular / polar / angle histogram (https://github.com/zifredder/CircHist), GitHub.
%   Adapted to a plotting function for ease of use. October 8, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make initial plot
fH = figure;
ax = polaraxes(fH);
obj1 = CircHist(sDist, nBins, 'parent', ax);
fH.Visible = 'on';

%% adjust appearance
% change color of bars
obj1.avgAngH.LineStyle = '--'; % make average-angle line dashed
obj1.avgAngH.LineWidth = 1; % make average-angle line thinner
obj1.colorAvgAng = [.5 .5 .5]; % change average-angle line color

% remove offset between bars and plot-center
rl = rlim(obj1.polarAxs); % get current limits
obj1.setRLim([0, rl(2)]); % set lower limit to 0

% draw circle at r == 0.5 (where r == 1 would be the outer plot edge)
rl = rlim(obj1.polarAxs);
obj1.drawCirc((rl(2) - rl(1)) /2, '--b', 'LineWidth', 2)
obj1.scaleBarSide = 'right'; % draw rho-axis on the right side of the plot
obj1.polarAxs.ThetaZeroLocation = 'right'; % rotate the plot to have 0° on the right side
obj1.setThetaLabel('Direction', 'bottomleft'); % add label

% draw resultant vector r as arrow
delete(obj1.rH)
obj1.drawArrow(obj1.avgAng, obj1.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')

% Change theta- and rho-axis ticks
obj1.polarAxs.ThetaAxis.MinorTickValues = []; % remove dotted tick-lines
% thetaticks(0:90:360); % change major ticks
% rticks(0:4:20); % change rho-axis tick-steps
obj1.drawScale; % update scale bar
end

