function tc_shuffle_dMan_homeRun(input, varargin)
%TC_SHUFFLE (dMan version)
%   As in Mimica et al. 2018
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

%% 0. Get stuff from dMan

% parse inputs
S = dManDe.Settings();
inp = inputParser();
inp.addParameter('goalLoc', 37)
% inp.addParameter('addPeakRate',false)
inp.parse(varargin{:});
inp.KeepUnmatched = true;
P = inp.Results;

% find the unit of interest?
% [unit, tracker] = dManDe.helpers.handleInput(input, varargin{:});
[unit, tracker] = dManDe.helpers.handleInput(input);

recType = unit.recording.sessions.type;   
if contains(recType, 'fm')

    % check for valid goal location
    if ~(P.goalLoc > 0 && P.goalLoc < 38)
        warning('first varargin (goal location) can only be a number between 1 and 37')
        P.goalLoc = 37;
    end  

    % grab position data
    t_trial = tracker.t; fs = mode(diff(t_trial));
    x_trial = tracker.x;
    y_trial = tracker.y;
    position_trial = [t_trial x_trial y_trial];
    hd_trial = rad2deg(tracker.hd); % 0 to 360 deg

    % grab goal location for this unit
    xGoal = unit.recording.sessions.userData.wellLocations(P.goalLoc,1);
    yGoal = unit.recording.sessions.userData.wellLocations(P.goalLoc,2);
    ref_point = [xGoal, yGoal];

    % get other information
    ST_unfiltered_trial = unit.t; % spike times (im assuming), in seconds

    % speed threshold the spike times (do we want this?)
    % lower_bound = 0.04;
    % s = unit.recording.trackers.speed; % grab speed (m/s)
    % ST_trial = speed_filter(t, s, ST_unfiltered_trial, lower_bound);

    % get FMEVNTS struct!!
    fmEvents = unit.recording.sessions.userData.fmEvents;

    % get homeRun struct (make sure this function is on path)
    [homeRun] = dManDe.helpers.getHomeRuns(position_trial, ST_unfiltered_trial, hd_trial, ref_point, fmEvents);

    % rename all the home events
    position = homeRun.position;
    hd = homeRun.hd;
    ST = homeRun.spiketimes;

    %% I. Calculate shuffled distribution

    % define the number of shuffles we want
    % we may need to change this based on # of spikes?
    total_shuffles = 100;

    % generate random floating-point values between a & b
    a = -60; b = 60;
    r = (b-a).*rand(total_shuffles*10,1) + a;

    % only keep values that are larger than 15(s) and smaller than -15(s)
    r_thresh = r(or(r>15, r<-15));

    % randomly sample from 'r_thresh' to get a random string of 'shift' values.
    shiftVals = datasample(r_thresh, total_shuffles);

    % clear variables
    tcVals_shift = [];

    num_tc_bins = 39; % number of tuning curve bins (this is hard-coded right now)

    for iter = 1:total_shuffles
        if ~isempty(ST) % check that there is at least one spike!
            % how much to shift (in s) for this iteration
            shift = shiftVals(iter);

            % get shifted timestamps
            ST_shift = dManDe.plots.circShift_TimeStamps(position, ST, shift);

            % get values for tuning curve (egoBear)
            [tcVals_shift(iter, :), ~] = egoBearing(position, ST_shift, hd, ref_point);
        else
            % if there are no spikes, make a tuning curve of NaNs
            tcVals_shift(iter, :) = ones(1,num_tc_bins)*NaN; 
        end
    end

    % take mean (along the columns)
    mean_tc = nanmean(tcVals_shift, 1);

    % take standard deviation (along the columns)
    std_tc = std(tcVals_shift, [], 1);

    % calculate the upper & lower bounds of y (+/- 2 stds)
    yu = mean_tc + 2*std_tc;
    yl = mean_tc - 2*std_tc;


    %% II. Calculate 'real' tuning curve
    [tcVals_real, binCtrs] = egoBearing(position, ST, hd, ref_point);


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
    title("Egocentric Bearing", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
    ylabel("fr (Hz)")
    xlim([0 360])
    xticks([0 90 180 270 360])
    xlabel("angle (deg)")
    box off

    hold off;
else
    disp('Open field session...')
end

end


