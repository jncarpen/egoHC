function ST_shftd = circShift_TimeStamps(position, ST, shiftVal)
%CIRCSHIFT_TIMESTAMPS Summary of this function goes here
%   Circularly shift timestamps (shiftVal is in seconds)
%   Inputs:
%   'position'          [t x y x2 y2]
%   'ST'                spiketimes array, with times in seconds
%   'shiftVal'          value to shift the spikes by (in seconds)
%   Outputs:
%   'ST_shftd'          circularly shifted timestamps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get timestamps
t = position(:,1);

tMax = nanmax(t); tMin = nanmin(t);

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

