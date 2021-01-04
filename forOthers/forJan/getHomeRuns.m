function [homeRun] = getHomeRuns(position, ST, hd, goal_loc, fmEvents)
%GETHOMERUNS Summary of this function goes here
%   Inputs:
%   'position'              Tx3 position vector in the form [t x y], where
%                           T is the number of position samples
%   'ST'                    Sx1 vector of spike times, where S is the
%                           number of spikes the neuron had in a given session.
%   'goal_loc'              1x2 vector with specified goal location
%   Outputs:
%   'HomeRuns'              struct.
%   J.Carpenter
%   Last modified: October 30, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse position
t = position(:,1); % time (s)
x = position(:,2);
y = position(:,3);
fs = mode(diff(t)); % sampling freq

% parse goal/reference location
xGoal = goal_loc(1,1);
yGoal = goal_loc(1,2);

% make spiketrain for cell (unsmoothed) *
% not sure if its less biased to smooth before or after?
startTime = t(1); stopTime = t(end);
ST = ST(ST < stopTime & ST > startTime);
edgesT = linspace(startTime,stopTime,numel(t)+1);
spiketrain = histcounts(ST,edgesT);

% grab info for this session (figure out how this is done in dMan)
eventStruct = fmEvents{1,sessNum}.fmEvents;
type = eventStruct.type;
times = eventStruct.times; % seconds?
locs = eventStruct.locations; % don't think we need this?

% grab home drinking events
type_count = count(type, "DRINKING_HOME");
home = times.*type_count; % pull out drinking_random events
home = home(any(home,2),:); % pulls out 0 values
home_ts = knnsearch(t, home(:,1));
home = t(knnsearch(t, home(:,1))); % match with closest timestamp


% remove first event IF its a 'drinking_home' event
if type{1,1} == 'DRINKING_HOME'
    home_ts = home_ts(2:end,:);
    home = home(2:end,:);
end

% grab decision point (only first column)
type_count = count(type, "DECISION_POINT");
decision = times.*type_count; % pull out drinking_random events
decision = decision(any(decision,2),:); % pulls out 0 values
decision_ts = knnsearch(t, decision(:,1)); 
decision = t(knnsearch(t, decision(:,1))); % match with closest timestamp


% check that decision + home vectors are same length
if length(home_ts) ~= length(decision_ts)
    disp('Error: home and decision vectors are of unequal length.')
    % add line where I fix the bug lol
end

% get position stamps
event_indices = [];
for evnt = 1:length(decision_ts)
    idx_now = decision_ts(evnt):1:home_ts(evnt);
    event_indices = [event_indices, idx_now];
end
event_indices = event_indices'; % transpose for your pleasure

event_timestamps = t(event_indices);
event_x = x(event_indices);
event_y = y(event_indices);
event_position = [event_timestamps, event_x, event_y];
event_hd = hd(event_indices);

% get spiketrain + smoothed spiketrain
event_spiketrain = spiketrain(event_indices); sigma = 2;
event_spiketrain_smooth = imgaussfilt(event_spiketrain, sigma, 'Padding', 'replicate');

% grab spiketimes that fall within events of interest
st_indices = [];
for evnt = 1:length(decision)
    st_indices = [st_indices; find(ST>decision(evnt) & ST<home(evnt))];
end
event_spiketimes = ST(st_indices);


%% put everything into a struct for later use
homeRun.indices = event_indices;
homeRun.position = event_position;
homeRun.hd = event_hd;
homeRun.spiketrain = event_spiketrain;
homeRun.spiketrain_smooth = event_spiketrain_smooth;
homeRun.spiketimes = event_spiketimes;

end

