function [taskPhz] = parseTask(events, pos)
    %TASKPHASE Bin timestamps by task phase for Foster Maze.
    %   Make sure that both events and timestamps are in SECONDS.

    %   INPUT:
    %   events struct:      col1 is time event began, col2 is time event ended
    %   pos:                [t x y x1 y2]
    
    
    %   Task phases defined as:
    %   (P1) Session onset
    %   (P2) Home drinking
    %   (P3) Foraging
    %   (P4) Random drinking
    %   (P5) Trip home
    %   (P6) Session offset
    
    
    %   OUTPUT:
    %   taskPhz:            matrix binned by task phase.
    
    %   Jordan Carpenter, 2020.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% (P2 and P4) Define drinking event times
    % define sampling frequency
    t = pos(:,1); % pull timestamps
    Fs = mode(diff(t)); % video sampling freq
    start = min(t);
    stop = max(t);
    
    % create binary vector of (P2) Home Drinking events
    type = count(events.type, "DRINKING_HOME");
    home = events.times.*type; % pull out drinking_random events
    home = home(any(home,2),:); % pulls out 0 values
    home = [t(knnsearch(t, home(:,1))), t(knnsearch(t, home(:,2)))]; % match with closest timestamp
    
    % create vector of (P4) Random drinking events
    rand = events.times.*~type; % pull out drinking_home events                                                      
    rand = rand(any(rand,2),:); % pulls out 0 values
    rand = [t(knnsearch(t, rand(:,1))), t(knnsearch(t, rand(:,2)))];


    %% (P1) taskOnset vector and (P6) taskOffset vector
    
    sessOnset = start:Fs:events.times(1,1)-Fs;
    ON_idx = knnsearch(t, sessOnset');
    sessOnset = t(ON_idx);
    
    sessOffset = events.times(end,2)+Fs:Fs:stop;
    OFF_idx = knnsearch(t, sessOffset');
    sessOffset = t(OFF_idx);
    
    
    %% (P3 & P5) Foraging and trip home events
    
    % Find start/end times for foraging & trip home events
    for eventNum = 1:length(events.type)
        if eventNum < length(events.type)-1
            if events.type{eventNum,1} == "DRINKING_HOME" && events.type{eventNum+1,1}== "DRINKING_RANDOM"
                Foraging(eventNum, 1:2) = [events.times(eventNum,2)+Fs, ; events.times(eventNum+1,1)-Fs];
            elseif events.type{eventNum,1} == "DRINKING_RANDOM" && events.type{eventNum+1,1}== "DRINKING_HOME"
                tripHome(eventNum, 1:2) = [events.times(eventNum,2)+Fs, ; events.times(eventNum+1,1)-Fs];
            else
                disp("Events repeated, task out of order.") % find a way to deal with this error.
            end  
        end
    end
   
    % pull out zero values in Foraging and tripHome vectors
    Foraging = Foraging(any(Foraging,2),:);
    Foraging = [t(knnsearch(t, Foraging(:,1))), t(knnsearch(t, Foraging(:,2)))]; % match with closest timestamp

    tripHome = tripHome(any(tripHome,2),:);
    tripHome = [t(knnsearch(t, tripHome(:,1))), t(knnsearch(t, tripHome(:,2)))]; % match with closest timestamp


    %% Put it all together...
    
    % bin event data to match timestamps
    edgesT = linspace(start,stop,numel(t)+1); % binsize is close to video frame rate (~30ms)
    
    homeTimes = [];
    for i = 1:length(home)
        homeList = home(i,1):Fs:home(i,2);
        homeTimes = [homeTimes, homeList];
    end
    HT_idx = knnsearch(t,homeTimes');
    homeTimes = t(HT_idx);
    
    randTimes = [];
    for i = 1:length(rand)
        randList = rand(i,1):Fs:rand(i,2);
        randTimes = [randTimes, randList];
    end
    R_idx = knnsearch(t,randTimes');
    randTimes = t(R_idx);

    foragingTimes = [];
    for i = 1:length(Foraging)
        foragingList = Foraging(i,1):Fs:Foraging(i,2);
        foragingTimes = [foragingTimes, foragingList];
    end
    F_idx = knnsearch(t,foragingTimes');
    foragingTimes = t(F_idx);
    
    tripHomeTimes = [];
    for i = 1:length(tripHome)
        tripHomeList = tripHome(i,1):Fs:tripHome(i,2);
        tripHomeTimes = [tripHomeTimes, tripHomeList];
    end
    TH_idx = knnsearch(t,tripHomeTimes');
    tripHomeTimes = t(TH_idx);
    
    %% Create a new vector- each timestamp has a number 1-6 corresponding
    % to task phase during that time.
    
    taskPhz = t; % copy time vector
    
    taskPhz(ON_idx) = 1; % session onset
    taskPhz(HT_idx) = 2; % drinking at home well
    taskPhz(F_idx) = 3; % foraging
    taskPhz(R_idx) = 4; % drinking at random well
    taskPhz(TH_idx) = 5; % trip home
    taskPhz(OFF_idx) = 6; % session offset

end

