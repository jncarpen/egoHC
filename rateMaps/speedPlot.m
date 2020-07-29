function speedPlot(SpikeTimes, pos, speed)
%SPEEDPLOT Summary of this function goes here
%   Rich (DMan)

t = pos(:,1); % pull timestamps

addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));
smap = analyses.map([t,speed], SpikeTimes,'binWidth', 0.05,'smooth', 2);
ax = gca;
plot(ax,smap.x(2:end),smap.z,'k', 'linewidth', 1.5);
ylim([0,ceil(max(smap.z)+1)]);
xlim([0,0.6]);
axis square
xlabel('speed (m/sec)');
ylabel('spike rate (Hz)');
box off
end

