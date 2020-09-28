function tc_stats_heatmap(pos_cm, hd, UniqueID, SpikeTimes_thresh, sessNum, unitNum)
    %TC_STATS_HEATMAP

    % get position vector (cm) now
    posNow = pos_cm{1,sessNum};
    timeNow = posNow(:,1);
    hdNow = hd{1,sessNum};
    STNow = SpikeTimes_thresh{1,sessNum}{1,unitNum};
    UID = UniqueID{1,sessNum}{1,unitNum};
    
    % pick a random normal distribution of spikes
%     STNow_sim = randi(floor(nanmax(timeNow)),length(STNow),1);
    STNow_sim = circshift(STNow,30); % shift by 30 seconds
    
    % set arbitrary values for binWidth and percentile
    binWidth = 18; % in degrees
    percentile = 95;

    % generate vector of reference locations (in cm)
    [refVec] = generate_reference_pnts(posNow);

    for refIdx = 1:length(refVec)
        
        % choose refLoc (and refLoc2)
        refLoc(1,1) = refVec(refIdx, 1); refLoc(1,2) = refVec(refIdx, 2);
        refLoc2(1,1) = refVec(refIdx, 1); refLoc2(1,2) = refVec(refIdx, 2);
                
        % generate egoBearing tuning curve- does this need to be in degrees
        % to work?
        [tcVals_egoAng] = egoBearing(pos_cm, STNow, refLoc, refLoc2, hd, sessNum, "False", "deg");
        [tcVals_egoAng_sim] = egoBearing(pos_cm, STNow_sim, refLoc, refLoc2, hd, sessNum, "False", "deg");

        % get tuning curve statistics
        tcStat = analyses.tcStatistics(tcVals_egoAng', binWidth, percentile);
        tcStat_sim = analyses.tcStatistics(tcVals_egoAng_sim', binWidth, percentile);

        
        % make matrix for resulting r values (mean vector length)
        r(refIdx) = tcStat.r;
        rSim(refIdx) = tcStat_sim.r;
%         peakDir(refIdx) = tcStat.peakDirection;
        
    end
    
    % reshape the matrix to be 10x10 and make into a heatmap
    start = 1;
    stop = 10;
    for ii = 1:10
        r_mat(:,ii) = flip(r(start:stop))';
        rSim_mat(:,ii) = flip(rSim(start:stop))';
%         peakDir_mat(:,ii) = flip(peakDir(start:stop))';
        start = start + 10; stop = stop + 10;
    end
   
    % make a heatmap (mean vector length- real data)
    figure
    imagesc(zscore(r_mat));
    figTit = strcat('MVL', ' ID', sprintf('%.f', UID), ' S', sprintf('%.f', sessNum));
    set(gca,'YDir','normal')
    pbaspect([1 1 1])
    title(figTit)
    colorbar
%     caxis([0, 1])
    
    % make a heatmap (mean vector length- sim data)
    figure
    imagesc(rSim_mat);
    figTit = strcat('MVLSIM', ' ID', sprintf('%.f', UID), ' S', sprintf('%.f', sessNum));
    set(gca,'YDir','normal')
    pbaspect([1 1 1])
    title(figTit)
    colorbar
    caxis([0, 1])
%     
%     % heatmap (peak direction)
%     imagesc(peakDir_mat);
%     set(gca,'YDir','normal')
%     pbaspect([1 1 1])
%     title("peakdir")
%     colorbar
%     colormap("hsv")
    
    
end
