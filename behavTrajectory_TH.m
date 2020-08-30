
function behavTrajectory_TH(events, pos)
    %BEHAVTRAJECTORY Generate plot of behavioral trajectory for trip home
    % events.
    %   Make sure that both events and timestamps are in SECONDS.
    %
    %   INPUT:
    %   events struct:              col1 is time event began, col2 is time event ended
    %   pos:                        [t x y x1 y2]
    %
    %
    %   Task phases defined as:
    %   (P1) Session onset
    %   (P2) Home drinking
    %   (P3) Foraging
    %   (P4) Random drinking
    %   (P5) Trip home
    %   (P6) Session offset
    %
    %
    %   OUTPUT:
    %   behavTrajectory plot:       script will generate a plot of all behavioral trajectories
    %                               from the time the animal *stops*
    %                               drinking from the random reward well
    %                               and takes its "trip home". This is
    %                               plotted in different colors for each
    %                               instance/trial.
    %   
    %
    %   Jordan Carpenter, 2020.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Define drinking event times
    
    % pull out important values
    t = pos(:,1); % pull timestamps
    x = pos(:,2);
    y = pos(:,3);
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
    
    
    %% Foraging and trip home events
    
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

    %% Plot behavioral trajectories
    
    % (1) Trip Home:
    figure
    for i = 1:length(tripHome)
        thPos = [];
        tripHomeList = tripHome(i,1):Fs:tripHome(i,2);
        TH_idx = knnsearch(t,tripHomeList');
        thPos(:,1) = x(TH_idx);
        thPos(:,2) = y(TH_idx);
        plot(thPos(:,1), thPos(:,2), 'LineWidth', 1.25)
        title("behavioral trajectory: trip home")
        xlabel("x (units)")
        ylabel("y (units)")
        hold on
        %pause
        %clf
    end
end
