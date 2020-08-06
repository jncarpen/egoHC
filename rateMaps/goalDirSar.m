function goalDirSar(pos, hwLoc, hd, SpikeTimes)
%SARELPLOT Summary of this function goes here
%   INPUTS
%   pos:            [t x y x2 y2]. 5xT vector, where T is number of
%                   timestamps.
%   hwLoc:          A single cell from the 1xS cell array called hwLoc, where 
%                   S is the number of sessions for a particular animal (hwLoc
%                   can be generated from the getWellLoc.m function). This will
%                   be a scalar value of the location of the home well for FM
%                   trials. Numbers will be in FM coordinates and *not* XY
%                   coordinates, for example: 36 or 37.
%   hd:             1xT vector of hd (in degrees) for the session, where T is the
%                   number of timestamps (should be same length as pos).
%
%   OUTPUTS
%
%
%   NOTE: For now, this is computed using HD, but once we compute *moving
%   direction*, it should be added to the function as an option.
%   
%
%   Jordan Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add BNT to path
    addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));
    
    % parse position vector
    t = pos(:,1);
    x = pos(:,2);
    y = pos(:,3);
    sampleRate = mode(diff(t));
   
    % manually define x/yGoal (for now)
    switch hwLoc
        case 36
            xGoal = 372.31;
            yGoal = 269.72;
        case 37
            xGoal = 372.58;
            yGoal = 305.75;            
    end
    
    % compute difference between currentPos and goalPos
    a = x - xGoal;
    b = y - yGoal;
    
    % calculate goal direction
    theta = atan2(b, a);
    goalDir = abs(theta - hd);
   
    % get spkGD
    idx = knnsearch(goalDir, SpikeTimes);
    spkGD = goalDir(idx);
    
    % calculate GD tuning curve
    tc_GD = analyses.turningCurve(spkGD, pos, sampleRate, 'binWidth', 10);
    
    % plot tuning curve
    plot(tc_GD(:,1), tc_GD(:,2), 'Color', 'k', 'LineWidth', 1.5)
    title("GoalDir (HD) TC")
    xlabel("angle to goal (deg)")
    ylabel("fr (Hz)")
    xlim([0 360])
    box off
    
end

