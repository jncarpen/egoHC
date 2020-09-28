function tc_stats_heatmap(pos_cm, hd, UniqueID, SpikeTimes_thresh, sessNum, unitNum)
    %TC_STATS_HEATMAP
    
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

    for refIdx = 1:length(refVec)

        % choose refLoc (and refLoc2); refLoc=refLoc2 for now
        refLoc(1,1) = refVec(refIdx, 1); refLoc(1,2) = refVec(refIdx, 2);
        refLoc2(1,1) = refVec(refIdx, 1); refLoc2(1,2) = refVec(refIdx, 2);

        % generate egoBearing tuning curve
        [tcVals_egoAng] = egoBearing(pos_cm, STNow, refLoc, refLoc2, hd, sessNum, "False", "deg");
        
        % get tuning curve statistics
        tcStat = analyses.tcStatistics(tcVals_egoAng', binWidth, percentile);

        % make matrix for resulting r values (mean vector length)
        r(refIdx) = tcStat.r;
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
    
    for iter = 1:numIter
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
            rSim(refIdx) = tcStat_sim.r;
            
            % calculate mean rSim value
            if refIdx == 1
                rSim_mean(refIdx) = rSim(refIdx);
            elseif refIdx > 1
                rSim_mean(refIdx) = (rSim_mean(refIdx) + rSim(refIdx))/2;
            end
            
        end

        % Reshape the matrix (10x10)
        start = 1; stop = 10;
        for ii = 1:10
            rSim_mat(:,ii) = flip(rSim(start:stop))';
            rSim_mean_mat(:,ii) = flip(rSim_mean(start:stop))';
            start = start + 10; stop = stop + 10;
        end
        
        % plot each iteration
        imagesc(rSim_mat-rSim_mean_mat)
        figTit = strcat('MVL', ' ID', sprintf('%.f', UID), ' S', sprintf('%.f', sessNum));
        set(gca,'YDir','normal')
        pbaspect([1 1 1])
        title(figTit)
        
        % add stuff to a matrix 
        r_sim_cell{1,iter} = rSim_mat;
        shiftVal = shiftVal + 2;
        
    end
    

    %% Plot
    
    % make a heatmap (mean vector length- real data)
    figure
    imagesc(zscore(r_mat-rSim_mat));
    figTit = strcat('MVL', ' ID', sprintf('%.f', UID), ' S', sprintf('%.f', sessNum));
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
