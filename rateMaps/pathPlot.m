function pathPlot(posx, posy, t_s, spikeTs)
%PATHPLOT: make simple pathplot

% make sure BNT is on the path
addpath(genpath('C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001'));

pos = [t_s, posx, posy];

% from BNT
[spkPos] = data.getSpikePositions(spikeTs', pos);


figure
plot(posx, posy, 'k')
hold on
scatter(spkPos(:,2),spkPos(:,3), 2, 'r')
title("Path Plot")
xlabel("X-Coordinate")
ylabel("Y-Coordinate")
return

end

