function [taskPhz] = parseTask(events, timestamps)
    %TASKPHASE Bin timestamps by task phase for Foster Maze.
    %   Make sure that both events and timestamps are in SECONDS.

    %   INPUT:
    %   events struct: col1 is time event began, col2 is time event ended

    %   Task phases defined as:
    %   (P1) Session onset
    %   (P2) Home drinking
    %   (P3) Foraging
    %   (P4) Random drinking
    %   (P5) Trip home
    %   (P6) Session offset
    
    %   Jordan Carpenter, 2020.
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% (P2 and P4) Define drinking event times
    % define sampling frequency
    sampleRate = mode(diff(timestamps));
    start = min(timestamps);
    stop = max(timestamps);
    
    % create binary vector of (P2) Home Drinking events
    type = count(events.type, "DRINKING_HOME");
    home = events.times.*type; % pull out drinking_random events
    home = home(any(home,2),:); % pulls out 0 values
    
    % create vector of (P4) Random drinking events
    rand = events.times.*~type; % pull out drinking_home events                                                      
    rand = rand(any(rand,2),:); % pulls out 0 values


    %% (P1) taskOnset vector and (P6) taskOffset vector
    
    sessOnset = start-0.2:sampleRate:events.times(1,1)-0.001;
    sessOffset = events.times(end,2)+0.001:sampleRate:stop;
    
    
    %% (P3 & P5) Foraging and trip home events
    
    % Find start/end times for foraging & trip home events
    for eventNum = 1:length(events.type)
        if eventNum < 57
            if events.type{eventNum,1} == "DRINKING_HOME" & events.type{eventNum+1,1}== "DRINKING_RANDOM"
                Foraging(eventNum, 1:2) = [events.times(eventNum,2)+0.001, ; events.times(eventNum+1,1)-0.02];
            elseif events.type{eventNum,1} == "DRINKING_RANDOM" & events.type{eventNum+1,1}== "DRINKING_HOME"
                tripHome(eventNum, 1:2) = [events.times(eventNum,2)+0.001, ; events.times(eventNum+1,1)-0.02];
            else
                disp("Events repeated, task out of order.") % find a way to deal with this error.
            end  
        end
    end
   
    % pull out zero values in Foraging and tripHome vectors
    Foraging = Foraging(any(Foraging,2),:);
    tripHome = tripHome(any(tripHome,2),:);

    %% Put it all together...
    
    % bin event data to match timestamps
    edgesT = linspace(start,stop,numel(timestamps)+1); % binsize is close to video frame rate (~30ms)
    
    homeTimes = [];
    for i = 1:length(home)
        homeList = home(i,1):sampleRate:home(i,2);
        homeTimes = [homeTimes, homeList];
    end
    
    randTimes = [];
    for i = 1:length(rand)
        randList = rand(i,1):sampleRate:rand(i,2);
        randTimes = [randTimes, randList];
    end

    foragingTimes = [];
    for i = 1:length(Foraging)
        foragingList = Foraging(i,1):sampleRate:Foraging(i,2);
        foragingTimes = [foragingTimes, foragingList];
    end
    
    tripHomeTimes = [];
    for i = 1:length(tripHomeTimes)
        tripHomeList = tripHome(i,1):sampleRate:tripHome(i,2);
        tripHomeTimes = [tripHomeTimes, tripHomeList];
    end
    
    % get rid of times outside of tracking times
    sessOnset = sessOnset(sessOnset < stop & sessOnset > start);
    homeTimes = homeTimes(homeTimes < stop & homeTimes > start);
    foragingTimes = foragingTimes(foragingTimes < stop & foragingTimes > start);
    randTimes = randTimes(randTimes < stop & randTimes > start);
    tripHomeTimes = tripHomeTimes(tripHomeTimes < stop & tripHomeTimes > start);
    sessOffset = sessOffset(sessOffset < stop & sessOffset > start);
   
    % bin everything
    [sOnCount,~,sOnInd]=      histcounts(sessOnset,edgesT);
    [homeCount,~,homeInd]=    histcounts(homeTimes,edgesT);
    [forCount,~,forInd]=      histcounts(foragingTimes,edgesT);
    [randCount,~,randInd]=    histcounts(randTimes,edgesT);
    [tripHCount,~,tripHInd]=  histcounts(tripHomeTimes,edgesT);
    [sOffCount,~,sOffInd]=    histcounts(sessOffset,edgesT);
    
    % Make sure there is no overlap 
    overlapVal = sum(find(sOnCount==1 & homeCount==1 & forCount==1 &
    randCount==1 & tripHCount==1 & sOffCount==1));
    
    % Squish into one vector
    % Here things are getting messed up (too many 0s and 9s)
    if overlapVal == 0
        taskPhz = (sOnCount + (2.*homeCount) + (3.*forCount) + (4.*randCount) + (5.*tripHCount) + (6.*sOffCount))';
    else
        disp("Error: Fix overlap in individual taskphase vectors.")
    end
    
    %% Misc. Error Messages
    
    if sum(taskPhz == 0) > 0
        disp("Error: Fix zeros in taskPhz vector.")
    elseif sum(taskPhz > 6) > 0
        disp("Error: Values over 6 present in taskPhz vector.")
    end
    
    return
end

