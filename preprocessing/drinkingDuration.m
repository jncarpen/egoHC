function [homeDur, randDur] = drinkingDuration(events)
%DRINKINGDURATION Compute duration of drinking periods
%   Jordan Carpenter, July 1, 2020.

duration = events.times(:,2) - events.times(:,1);
type = count(events.type, "DRINKING_HOME");

homeDur = duration.*type;
homeDur = homeDur(any(homeDur,2),:); % pulls out 0 values

randDur = duration.*~type;                                                            
randDur = randDur(any(randDur,2),:); % pulls out 0 values
