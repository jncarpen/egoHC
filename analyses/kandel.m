%% KANDEL ANALYSIS
%  October 2020.

%% I. Get data
% A. For Jan Sigurd's data:

% choose a session/unit number
sessNum = 27; % unitNum = 4;

% get position data for session
position = pos_cm{1,sessNum};

% choose reference points
ref_point = hwCoord{1,sessNum};
ref_point2 = boxCtr{1,sessNum};

% get bin center locations
nBins = 10; % divide the arena into 100 2D spatial bins
[binCtrs] = get_spatial_bin_centers(position, nBins);

% define 'angle_of_interest'
angle_of_interest = 90;

% simulate an egocentric bearing cell with behavioral data
[S, SpkTrn, head_direction] = simulate_ego_cell(position, ref_point, angle_of_interest);
SpkTrn = SpkTrn';

trialType{1,sessNum}


% % (for real data):
% head_direction = get_hd(position);
% [spkTrainNow_thresh] = spkTrain_thresh(pos_cm, SpikeTimes_thresh, SpikeTrain); % get speed-thresholded spiketrain
% SpkTrn = spkTrainNow_thresh{1,sessNum}{1,unitNum}';
% S = SpikeTimes_thresh{1,sessNum}{1,unitNum}; % spikeTimes (thresholded)

%% II. Generate spike plots

figure
set(gcf,'color','w');
hold on;
pathPlot_hd(position, S, head_direction)
h1 = plot(ref_point(1,1), ref_point(1,2), 'o', 'MarkerSize', 12);
set(h1, 'markerfacecolor', 'k');
% h2 = plot(ref_point2(1,1), ref_point2(1,2), 'o', 'MarkerSize', 12);
% set(h2, 'markerfacecolor', 'red');
% legend("path", "spikes", "other", "goal loc", "Location", "northwestoutside")
hold off;

%% III. Calculate (+plot) tuning curves + statistics

% make binned occupancy plots

[hd_occ, allo_occ, ego_occ, time_occ] = get_binned_occupancy(position, ref_point, "ego");
sgtitle("Occupancy(Egocentric Bearing): Session 27")


% run 'TC_stats_2DBins.m' for a particular reference point (refPnt)
[HD_TC, ALLO_TC, EGO_TC, HD_ST, ALLO_ST, EGO_ST, allo_angleAtPeak, ego_angleAtPeak, allo_circVar, ego_circVar, ctrLocs] = TC_stats_2DBins(position, head_direction, SpkTrn, ref_point, 1);

% get scores for this reference points
[HD_mean_stats, HD_sum_stats] = score_tuning_curve(HD_ST);
[ALLO_mean_stats, ALLO_sum_stats] = score_tuning_curve(ALLO_ST);
[EGO_mean_stats, EGO_sum_stats] = score_tuning_curve(EGO_ST);

% mean vector length
ego_test_MVL = EGO_sum_stats.MVL;
allo_test_MVL = ALLO_sum_stats.MVL;
hd_test_MVL = HD_sum_stats.MVL;

% mean direction
ego_test_MD = EGO_mean_stats.mean_direction;
allo_test_MD = ALLO_mean_stats.mean_direction;
hd_test_MD = HD_mean_stats.mean_direction;

% peak rate
ego_test_pr = EGO_sum_stats.peak_rate;
allo_test_pr = ALLO_sum_stats.peak_rate;
hd_test_pr = HD_sum_stats.peak_rate;

fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(fig,'color','w');
set(fig, 'PaperPositionMode', 'auto')
sgtitle("TITLE")
cellNum = 1;
for row = 1:10
    for col = 1:10
        map_axis = deg2rad(EGO_TC{row,col}(:,1));
        tc_vals = EGO_TC{row,col}(:,2);
        
        % get mean vector length
        binWidth = 9; % should be same as that in 'TC_stats_2DBins'
        percentile = 95; % this is arbitary for the information we're looking at
        tcStat = analyses.tcStatistics(EGO_TC{row,col}, binWidth, percentile);
        tcStat_hd = analyses.tcStatistics(HD_TC{row,col}, binWidth, percentile);

        % make scale vector (for quiver plot)
        MVL_ego(cellNum) = tcStat.r;
        MVL_hd(cellNum) = tcStat_hd.r;
        peakDir_ego(cellNum) = tcStat.peakDirection;
        peakDir_hd(cellNum) = tcStat_hd.peakDirection;
        mean_ego(cellNum) = tcStat.mean; % mean direction
        mean_hd(cellNum) = tcStat_hd.mean;

        % name each plot
%         plot_title = strcat('MVL:', sprintf('%.4f', tcStat.r), '\PR:', sprintf('%.2f', nanmax(tc_vals)));
        plot_title = strcat('MVL:', sprintf('%.4f', tcStat.r), '\PD:', sprintf('%.f', peakDir_ego(cellNum)));
        
        % make it wrap around (?? why did I do this again?)
        map_axis = [map_axis;0]; 
        tc_vals = [tc_vals; tc_vals(1)];
        
%         % plot tuning curves
%         subplot(10,10,cellNum) 
%         polarplot(map_axis,tc_vals, 'LineWidth', 1.10, 'color', 'red')
%         ax = gca;
%         ax.ThetaTick = [0 90 180 270];
%         ax.RTick = [];
%         title(plot_title)

        % Occupancy plots
        
        

        cellNum = cellNum + 1;
    end
end

%% IV. Generate quiver plots

% define 'theta' and 'scale' inputs
theta = mean_ego; scale = MVL_ego;
theta2 = mean_hd; scale2 = MVL_hd;
point_locations = binCtrs;
refPnt = ref_point;

% plot quiver (egoBearing)
[quiver_egoBearing] = plot_quiver_scaled(theta, point_locations, scale, refPnt);
title("Egocentric Bearing")

% plot quiver (head direction)
[quiver_hd] = plot_quiver_scaled(theta2, point_locations, scale2, refPnt);
title("Head Direction")
