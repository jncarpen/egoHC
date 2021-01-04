function  pathPlot_foraging(input, varargin)
%PATHPLOT 

[unit, tracker, ~, P, ~] = dManDe.helpers.handleInput(input, varargin{:});

recType = unit.recording.sessions.type;   
if contains(recType, 'fm')

    if P.speedFilter
        [vUnit, ~, ~] = tracker.speedFilter(unit.t,[P.minSpeed P.maxSpeed],'between');
        [vTracker, ~, ~] = tracker.speedFilter(tracker.t,[P.minSpeed P.maxSpeed],'between');
    else
        vUnit = true(numel(unit.t),1);
        vTracker = true(numel(tracker.t),1);
    end


    % grab unit information
    t_trial = tracker.t(vTracker);
    x_trial = tracker.x(vTracker);
    y_trial = tracker.y(vTracker);
    hd_trial = tracker.hd(vTracker);
    position_trial = [t_trial x_trial y_trial];
    md = tracker.movingDirection(vTracker);
    ST_unfiltered_trial = unit.t(vUnit); % raw spiketimes

    % grab goal location for this unit
    xGoal = unit.recording.sessions.userData.wellLocations(P.goalLoc,1);
    yGoal = unit.recording.sessions.userData.wellLocations(P.goalLoc,2);
    ref_point = [xGoal, yGoal];
    
    % get FMEVNTS struct!!
    fmEvents = unit.recording.sessions.userData.fmEvents;
    
    % get homeRun struct (make sure this function is on path)
    [foraging] = dManDe.helpers.getForaging(position_trial, ST_unfiltered_trial, hd_trial, ref_point, fmEvents);

    % rename all the home events
    position = foraging.position;
    t = position(:,1);
    x = position(:,2);
    y = position(:,3);
    hd = foraging.hd;
    tSpk = foraging.spiketimes;
    

    % speed threshold the spike times
    lower_bound = 0.00;
    s_trial = unit.recording.trackers.speed; % grab speed (m/s)
    s = s_trial(foraging.indices);
    tSpk = speed_filter(t, s, tSpk, lower_bound);

    if isempty(P.parent)
        P.parent = gca;
    end


    % Restrict spike angles if necessary
    if ~isempty(P.spikeAngleRange)
        if ismac
            spikeIdx = knnsearch(tSpk, t);
        else
            spikeIdx = binarySearch(tSpk, t);
        end

        if P.spikeAngleRange(2) < P.spikeAngleRange(1)
            P.spikeAngleRange(2) = P.spikeAngleRange(2) + 2*pi;
        end

        spikeAngles = hd(spikeIdx);
        spikeAngles = mod(spikeAngles - P.spikeAngleRange(1), 2*pi) + P.spikeAngleRange(1);
        spikesVAngle = spikeAngles >= P.spikeAngleRange(1) & spikeAngles <= P.spikeAngleRange(2);
        tSpk = tSpk(spikesVAngle);
    end

    % Interpolate spike positions (should produce cleaner results than finding
    % position at nearest pos sample)
    spikePosX = interp1(t, x, tSpk);
    spikePosY = interp1(t, y, tSpk);

    % Add small amount of jitter to prevent exact position overlap
    jittX = diff(prctile(x, [5 95])) * 0.001;
    jittY = diff(prctile(x, [5 95])) * 0.001;
    spikePosX = spikePosX + (rand(size(spikePosX)) - 0.5)*(jittX);
    spikePosY = spikePosY + (rand(size(spikePosY)) - 0.5)*(jittY);

    % Auto-calculate marker size based on spike number
    if strcmpi(P.markerSize, 'auto') || strcmpi(P.pathPlotLineWidth, 'auto')

        % Calculate axes area in points
        units0 = P.parent.Units;
        P.parent.Units = 'points';
        axArea = prod(P.parent.Position([3 4]));
        P.parent.Units = units0;

        if strncmpi(P.markerSize, 'auto', 4)

            % Parse optional user-specified marker size multiplier
            % (e.g. "auto3", "auto0.5")
            if length(P.markerSize) > 4
                multiplier = str2double(P.markerSize(5:end));
            else
                multiplier = 1;
            end

            % The desired spike area is proportional to the axes area
            nSpikes = length(tSpk);
            spikeArea = axArea / nSpikes * multiplier;

            % If plotting with Line, must convert spike area to radius
            % Scale factor of 12 converts the Scatter 'o' marker scale to the
            % Line '.' marker
            if useScatter
                P.markerSize = spikeArea;
            else
                P.markerSize = sqrt(spikeArea)/pi * 12;
            end
        end

        if strcmpi(P.pathWidth, 'auto')
            P.pathWidth = sqrt(axArea) / 250;
        end
    end

    if P.showPath
        for trial = 1:length(foraging.indices_split)
            h.path = line( ...
                x_trial(foraging.indices_split{1,trial}), y_trial(foraging.indices_split{1,trial}), ...    
                'parent', P.parent, ...    
                'color', P.pathPlotColour, ...
                'lineWidth', P.pathPlotLineWidth);
            hold on;
        end
    end
    hold(P.parent, 'on')
    if P.showBox
        if tracker.transformToTemplateBox
            boxCorners = tracker.boxTemplateCorners;
        else
            boxCorners = tracker.boxCorners;
        end
        x = boxCorners([1:end 1], 1);
        y = boxCorners([1:end 1], 2);

        % Get box centre of mass and apply scaling factor
        boxCom = mean(boxCorners);
        x = x + (x - boxCom(1))*(P.boxScalingFactor - 1);
        y = y + (y - boxCom(2))*(P.boxScalingFactor - 1);
        line(x, y, 'parent', P.parent, 'color', P.boxColor, 'lineWidth', P.boxLineWidth);
    end


    if P.useScatter 
        useColMapping = true;

        if strcmpi(P.spikeColors, 'thetaPhase')
            lfp = unit.electrodeData.loadLfpData(P.session);
            lfp.resample(100).filterBandPass([5 10]);
            idx = dMan.helpers.binarySearch(tSpk, lfp.t);
            colVar = mod(angle(hilbert(lfp.y(idx))), 2*pi);
        elseif any(strcmpi(P.spikeColors, {'md', 'movingDirection'}))
            idx = dMan.helpers.binarySearch(tSpk, t);
            colVar = md(idx);  
        elseif any(strcmpi(P.spikeColors, {'hd', 'headDirection'}))
            idx = dMan.helpers.binarySearch(tSpk, t);
            colVar = hd(idx);  
        else
            useColMapping = false;
        end

        if useColMapping
            spikeCol = mapcolors(colVar, [0 2*pi], 'hsv');
        else
            spikeCol = P.spikeColors;
        end

        if P.markerFilled
            extraArgs = {'filled'};
        else
            extraArgs = {};
        end

        
        h.spikes = scatter(P.parent, spikePosX, spikePosY, P.markerSize, spikeCol, extraArgs{:});

    else
        h.spikes = line(spikePosX, spikePosY, ...
            'lineStyle', 'none', ...
             'marker', '.', ...
             'color', P.pathPlotSpikeColour, ...
             'markerSize', P.pathPlotMarkerSize*1.5, ...
             'parent', P.parent);
    end

    axis(P.parent, 'equal', 'off', 'tight')
else
    disp('Open field session...')
end

end

% colorbar(hsv)

