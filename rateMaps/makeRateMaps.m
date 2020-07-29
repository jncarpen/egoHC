% MAKE RATEMAPS


% WHAT TO MAKE:
% Firing rate map (w/ peak FR)
% Spatial occupancy
% Speed occupancy
% acceleration occupancy

for i = 1:length(SpikeTimes{1,28})
map = analyses.map(pos{1,28}, SpikeTimes{1,28}{1,i}, 'smooth',15, 'binWidth',2.5);subplot(2,2,1)
plot.colorMap(map.z)
colormap(jet)
subplot(2,2,2)
pathPlot(pos{1,28},SpikeTimes{1,28}{1,i})
pause
clf
end

% saveLoc = "D:\Data\rateMaps";