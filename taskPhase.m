function [tpMat] = taskPhase(events, timestamps)
%TASKPHASE Bin timestamps by task phase.
%   Make sure that both events and timestamps are in SECONDS.

%   INPUT:
%   events struct: col1 is time event began, col2 is time event ended

%   Task phases defined as:
%   (1) Home drinking
%   (2) Foraging
%   (3) Random drinking
%   (4) Trip home

% define sampling frequency
sampleRate = mode(diff(timestamps));
start = min(timestamps);
stop = max(timestamps);

type = count(events.type, "DRINKING_HOME");

home = events.times.*type;
home = home(any(home,2),:); % pulls out 0 values

% list of all times when animal was drinking from home well
homeTimes = [];
for i = 1:length(home)
    metab = home(i,1):sampleRate:home(i,2);
    homeTimes = [homeTimes, metab];
end

rand = events.times.*~type;                                                            
rand = rand(any(rand,2),:); % pulls out 0 values

% list of all times when animal was drinking from random well
randTimes = [];
for i = 1:length(rand)
    metab = rand(i,1):sampleRate:rand(i,2);
    randTimes = [randTimes, metab];
end

% We now need to parse out times for the other two task phases \\

end

