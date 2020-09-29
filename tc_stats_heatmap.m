function tc_stats_heatmap(pos_cm, hd, UniqueID, SpikeTimes_thresh, sessNum, unitNum)
    %TC_STATS_HEATMAP
    %   INPUTS
    %   'pos_cm'                position of animal in centimeters, in the
    %                           form of a 1xS cell array, where S is the number
    %                           of sessions. Each cell is populated with a Tx5 matrix,
    %                           where T is the number of timestamps- looks like this:
    %                           [t x1 y2 x2 y2].
    %
    %   'hd'                    cell array of *head direction values* for a single
    %                           animal. Same format as 'pos_cm'.
    %
    %   'UniqueID'              cell array of UniqueID values for the animal.
    %   
    %   'SpikeTimes_thresh'     cell array of thresholded spike times (only
    %                           times when the animal was moving >5cm/s are
    %                           included.
    %
    %   'sessNum'               session number that you're interested in
    %                           looking at.
    %
    %   'unitNum'               unit number (not uniqueID though)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% grab information + set values
    
    posNow = pos_cm{1,sessNum};
    timeNow = posNow(:,1);
    hdNow = hd{1,sessNum};
    STNow = SpikeTimes_thresh{1,sessNum}{1,unitNum};
    UID = UniqueID{1,sessNum}{1,unitNum};
    binWidth = 18; % deg
    percentile = 95;

    
    %% Compute statistics for *real* data
    
    numIter = 30;
    [refVec] = generate_reference_pnts(posNow);
    r = zeros(length(refVec),1);
    
    for refIdx = 1:length(refVec)

        % choose refLoc (and refLoc2); refLoc=refLoc2 for now
        refLoc(1,1) = refVec(refIdx, 1); refLoc(1,2) = refVec(refIdx, 2);
        refLoc2(1,1) = refVec(refIdx, 1); refLoc2(1,2) = refVec(refIdx, 2);

        % generate egoBearing tuning curve
        [tcVals_egoAng] = egoBearing(pos_cm, STNow, refLoc, refLoc2, hd, sessNum, "False", "deg");
        
        % get tuning curve statistics
        tcStat = analyses.tcStatistics(tcVals_egoAng', binWidth, percentile);

        % make matrix for resulting r values (mean vector length)
        r(refIdx,1) = tcStat.r;
    end

    % Reshape the matrix (10x10)
    start = 1; stop = 10;
    for ii = 1:10
        r_mat(:,ii) = flip(r(start:stop))';
        start = start + 10; stop = stop + 10;
    end
    
    
    %% Compute statistics for *shuffled* data
    
    shiftVal = 30; % s
    r_sim_cell = cell(1,numIter);
    rSim = zeros(length(refVec),1); rSim_mean = zeros(length(refVec),1);
    
    for iter = 1:2%:numIter
        STNow_sim = circShift_TimeStamps(pos_cm, SpikeTimes_thresh, sessNum, unitNum, shiftVal);
        
        for refIdx = 1:length(refVec)

            % choose refLoc (and refLoc2); refLoc=refLoc2 for now
            refLoc(1,1) = refVec(refIdx, 1); refLoc(1,2) = refVec(refIdx, 2);
            refLoc2(1,1) = refVec(refIdx, 1); refLoc2(1,2) = refVec(refIdx, 2);

            % generate egoBearing tuning curve
            [tcVals_egoAng_sim] = egoBearing(pos_cm, STNow_sim, refLoc, refLoc2, hd, sessNum, "False", "deg");

            % get tuning curve statistics
            tcStat_sim = analyses.tcStatistics(tcVals_egoAng_sim', binWidth, percentile);


            % make matrix for resulting r values (mean vector length)
            rSim(refIdx,1) = tcStat_sim.r;
            
            % calculate mean rSim value
%             if refIdx == 1
%                 rSim_mean(refIdx,1) = rSim(refIdx,1);
%             elseif refIdx > 1
%                 rSim_mean(refIdx,1) = (rSim_mean(refIdx,1) + rSim(refIdx,1))/2;
%             end
            
        end

        % Reshape the matrix (10x10)
        start = 1; stop = 10;
        for ii = 1:10
            rSim_mat(:,ii) = flip(rSim(start:stop))';
%             rSim_mean_mat(:,ii) = flip(rSim_mean(start:stop))';
            start = start + 10; stop = stop + 10;
        end
        
        % add stuff to a matrix 
%         r_sim_cell{1,iter} = rSim_mat;
        r_sim_3D(:,:,iter) = rSim_mat;
        shiftVal = shiftVal + 2;
    end
    
    % take mean of the r_sim_3D matrix in the 3rd dimension (element-wise)
    meanMatrix = mean(r_sim_3D,3);

    %% Plot
    
    % make a heatmap (mean vector length- real data)
    figure
    imagesc(r_mat-meanMatrix);
    figTit = strcat('MVL(real)-mean(MVL(shuff))', ' ID', sprintf('%.f', UID), ' S', sprintf('%.f', sessNum));
    set(gca,'YDir','normal')
    pbaspect([1 1 1])
    title(figTit)
    colorbar
    
end


%% SCRATCH
%     Alternative ways to shift spikes:
%
%     (1)
%     pick a random normal distribution of spikes
%     STNow_sim = randi(floor(nanmax(timeNow)),length(STNow),1);
%
%     (2)
%     circularly shift spikes by 150 *bins*
%     STNow_sim = circshift(STNow,150); % shift by shiftVal bins
