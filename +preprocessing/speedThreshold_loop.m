function [SpikeTimes_thresh] = speedThreshold_loop(pos_cm, speed_cm, SpikeTimes)
%SPEEDTHRESHOLD_LOOP Summary of this function goes here
%   Detailed explanation goes here

SpikeTimes_thresh = cell(1,length(pos_cm));

for sessNum = 1:length(pos_cm)
    sessionArray = cell(1,length(SpikeTimes{1,sessNum}));
    
    % grab stuff for current session
    speedNow = speed_cm{1,sessNum};
    posNow = pos_cm{1,sessNum};
    
    % grab speed for first LED
    spdLED1 = speedNow(:,1);
    minSpd = nanmin(spdLED1);
    maxSpd = nanmax(spdLED1);
    t = posNow(:,1);
    sampleRate = mode(diff(t));
    
    for unitNum = 1:length(SpikeTimes{1,sessNum})
        ST_now = SpikeTimes{1,sessNum}{1,unitNum};

        % find closest timestamp & store values
        idx = knnsearch(t, ST_now);
        spkSpd = spdLED1(idx);
        greater_than_5_idx = find(spkSpd>5);
        speed_thresh_spikeTimes = ST_now(greater_than_5_idx);

        sessionArray{1,unitNum} = speed_thresh_spikeTimes;
    end
    SpikeTimes_thresh{1,sessNum} = sessionArray;
end


end

