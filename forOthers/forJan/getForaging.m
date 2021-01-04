function [foraging] = getForaging(position, ST, hd, goal_loc, fmEvents)
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
type_raw = fmEvents.type;
times = fmEvents.times; % seconds?


event_type = type_raw;
event_type_2 = type_raw;
    for i = 1:length(type_raw)-1
        if isequal(type_raw{i,1}, type_raw{i+1,1})
            event_type{i+1,1} = [];
            event_type_2{i+1,1} = 'REPEATED_EVENT';
        end  
    end

% get rid of empty cells
index = cellfun(@isempty, event_type) == 0;
event_type = event_type(index);

% grab home drinking events
type_count = count(event_type_2, "DRINKING_HOME");
home_drink = times.*type_count; % pull out drinking_random events
home_drink = home_drink(any(home_drink,2),:); % pulls out 0 values
home_ts = knnsearch(t, home_drink(:,2)); % grab end (col2) of homing event
home_drink = t(knnsearch(t, home_drink(:,2))); % match with closest timestamp


% grab *onset* of random drinking event (col1)
type_count = count(event_type_2, "RANDOM_DRINKING");
rand_drink = times.*type_count; % pull out drinking_random events
rand_drink = rand_drink(any(rand_drink,2),:); % pulls out 0 values
rand_ts = knnsearch(t, rand_drink(:,1)); 
rand_drink = t(knnsearch(t, rand_drink(:,1))); % match with closest timestamp


% remove first event IF its a 'drinking_random' event (this might cause
% some issues)
if event_type{1,1} == "DRINKING_RANDOM"
    rand_ts = rand_ts(2:end,:);
    rand_drink = rand_drink(2:end,:);
else
end

% check that decision + home vectors are same length
if length(home_ts) ~= length(rand_ts)
    disp('Error: home and random foraging vectors are of unequal length.')
    % add line where I fix the bug lol
end

% get position stamps
event_indices = [];
event_indices_cell = cell(1, length(rand_ts));
for evnt = 1:length(rand_ts)
    idx_now = home_ts(evnt):1:rand_ts(evnt);
    event_indices = [event_indices, idx_now];
    event_indices_cell{1, evnt} = idx_now;
end
event_indices = event_indices'; % transpose for your pleasure :D

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
for evnt = 1:length(rand_drink)
    st_indices = [st_indices; find(ST>home_drink(evnt) & ST<rand_drink(evnt))];
end
event_spiketimes = ST(st_indices);


%% put everything into a struct for later use
foraging.indices = event_indices;
foraging.indices_split = event_indices_cell;
foraging.position = event_position;
foraging.hd = event_hd;
foraging.spiketrain = event_spiketrain;
foraging.spiketrain_smooth = event_spiketrain_smooth;
foraging.spiketimes = event_spiketimes;

end

