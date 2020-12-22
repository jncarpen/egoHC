function tc_shuffle(position, ST, ref_point)
%TC_SHUFFLE 
%   As in Mimica et al. 2018
%   Inputs:
%   'position'      tx5 array, where t is the number of timestamps in the
%                   session. It should be in the form [t x y x2 y2].
%   'ST'            sx1 array, where s is the number of spikes the cell
%                   fired throughout the session. spike time values should
%                   be in *seconds*.
%   'ref_point'     reference point for calculation of egocentric bearing;
%                   should be in the form [ref_x, ref_y]
%   Output:
%   This function will output the plot of the egocentric bearing tuning
%   curve with the shuffled distribution.
%
% Alg:
% I shift the data n times (circularly, between +/-15-60 seconds) and make 
% n+1 tuning curves (of egocentric bearing). The +1 is the tuning curve of 
% the actual data, which I plot in red. Then I take the mean (at each point/bin)
% of the n shuffled tuning curves and plot, in grey, the mean +/- 2 standard deviations. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the example...
% ref_point = [89, 83];

%% I. Calculate shuffled distribution

% define the number of shuffles we want
total_shuffles = 100;

% tell the 'plot_egoBearing.m' function not to plot
doPlot = "False";
deg_or_rad = "deg";

% generate random floating-point values between a & b
a = -60; b = 60;
r = (b-a).*rand(total_shuffles*5,1) + a;

% only keep values that are larger than 15(s) and smaller than -15(s)
r_thresh = r(or(r>15, r<-15));

% randomly sample from 'r_thresh' to get a random string of 'shift' values.
shiftVals = datasample(r_thresh, total_shuffles);

% clear variables
tcVals_shift = [];

nBins = 39;

for iter = 1:total_shuffles
    % how much to shift (in s) for this iteration
    shift = shiftVals(iter);
    
    % get shifted timestamps
    ST_shift = circShift_TimeStamps(position, ST, shift);
    
    % get values for tuning curve (egoBear)
%     [tcVals_shift(iter,:), ~] = plot_egoBearing(position, ST_shift, ref_point, doPlot); % speed-thresholded (5cm/s)
    [tcVals_shift(iter,:), ~] = egoBearing(position, ST_shift, ref_point, ref_point, nBins, doPlot, "deg");
end

% take mean (along the columns)
mean_tc = nanmean(tcVals_shift, 1);

% take standard deviation (along the columns)
std_tc = std(tcVals_shift, [], 1);

% calculate the upper & lower bounds of y (+/- 2 stds)
yu = mean_tc + 2*std_tc;
yl = mean_tc - 2*std_tc;


%% II. Calculate 'real' tuning curve
[tcVals_real, binCtrs] = plot_egoBearing(position, ST, ref_point, doPlot); % speed-thresholded (5cm/s)
binCtrs

%% III. Plot tuning curves
% figure
% set(gcf,'color','w');

% plot shuffled data
fill([binCtrs fliplr(binCtrs)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none')
hold all
plot(binCtrs, mean_tc, ':k')

% plot actual data
plot(binCtrs, tcVals_real, 'Color', 'r', 'LineWidth', 1.5)

% format the plot
% title("Egocentric Bearing")
ylabel("fr (Hz)", 'FontSize', 20)
pbaspect([1 1 1])
xlim([0 360])
xticks([0 90 180 270 360])
xlabel("angle (deg)")
set(gca,'FontSize',20, 'FontName', 'Calibri Light')
box off

hold off;
end

