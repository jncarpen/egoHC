function [fig, map] = plot_ratemap(P, ST)
%PLOT_RATEMAP using BNT stuff

fig = figure;
map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 150/10); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
colormap(gca,'jet')
c2 = colorbar; c2.FontSize = 25;
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
box off

end

