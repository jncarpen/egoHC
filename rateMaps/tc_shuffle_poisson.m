function tc_shuffle_poisson(P, ST, RP, size_simdataset)
%TC_SHUFFLE 
%   As in Mimica et al. 2018
%   Inputs:
%   'position'          tx5 array, where t is the number of timestamps in the
%                       session. It should be in the form [t x y x2 y2].
%   'ST'                sx1 array, where s is the number of spikes the cell
%                       fired throughout the session. spike time values should
%                       be in *seconds*.
%   'ref_point'         reference point for calculation of egocentric bearing;
%                       should be in the form [ref_x, ref_y]
%
%   'total_shuffles'    scalar value of how big you want the simulated dataset
%                       to be.
%   Output:
%   This function will output an egocentric bearing tuning curve (of the
%   real data) in red and the mean tuning curve of a distribution of simulated cells
%   (using Seb's wAnG method) +/- 2 standard deviations, shown in grey
%   shading.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I. Calculate shuffled distribution


% tell the 'plot_egoBearing.m' function not to plot
doPlot = "False";

% break down 'refLoc'
rlX = RP(1,1);
rlY = RP(1,2);

% break down position data
x = P(:,2); y = P(:,3);
t = P(:,1); hd = get_hd(P);

% get spike positions 
clear spikes
spikes = zeros(length(ST), 3);
spk_idx  = knnsearch(t, ST);
spikes(:,1) = ST; spikes(:,2) = x(spk_idx);
spikes(:,3) = y(spk_idx);

for iter = 1:size_simdataset
    [simulatedSpikes,~] = spatialPoissonSimulation(P, spikes);
    [tcVals_shift(iter,:), ~] = plot_egoBearing(P, simulatedSpikes(:,1), RP, doPlot);
end

% take mean (along the columns)
mean_tc = nanmean(tcVals_shift, 1);

% take standard deviation (along the columns)
std_tc = std(tcVals_shift, [], 1);

% calculate the upper & lower bounds of y (+/- 2 stds)
yu = mean_tc + 2*std_tc;
yl = mean_tc - 2*std_tc;


%% II. Calculate 'real' tuning curve
 [tcVals_data, binCtrs] = plot_egoBearing(P, ST, RP, doPlot);


%% III. Plot tuning curves
% figure
% set(gcf,'color','w');

% plot shuffled data
fill([binCtrs fliplr(binCtrs)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none')
hold all
plot(binCtrs, mean_tc, ':k')

% plot actual data
plot(binCtrs, tcVals_data, 'Color', 'r', 'LineWidth', 1.5)

% format the plot
% title("Egocentric Bearing")
ylabel("fr (Hz)")
pbaspect([1 1 1])
xlim([0 360])
xticks([0 90 180 270 360])
xlabel("angle (deg)")
set(gca,'FontSize',20, 'FontName', 'Calibri Light')
box off

hold off;
end

