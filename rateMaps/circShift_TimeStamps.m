function ST_shftd = circShift_TimeStamps(pos_cm, SpikeTimes_thresh, sessNum, unitNum, shiftVal)
%CIRCSHIFT_TIMESTAMPS Summary of this function goes here
%   Circularly shift timestamps (shiftVal is in seconds)

t = pos_cm{1,sessNum}(:,1);
tMax = nanmax(t); tMin = nanmin(t);
ST = SpikeTimes_thresh{1,sessNum}{1,unitNum};

ST_shftd = [];
for spk = 1:length(ST)
    testShift = ST(spk) + shiftVal;
    if testShift > tMax
        diff = testShift - tMax;
        ST_shftd = [ST_shftd, tMin+diff];
    else
        ST_shftd = [ST_shftd, testShift];
    end
end

% sort values in ascending order
ST_shftd = sort(ST_shftd', 'ascend');
end

