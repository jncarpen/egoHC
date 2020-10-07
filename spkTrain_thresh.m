function [SpikeTrain_thresh] = spkTrain_thresh(position, SpikeTimes, SpikeTrain)
%SPIKETRAIN_THRESHOLDED Summary of this function goes here
%   INPUTS
%   'position'      position in centimeters
%   'SpikeTimes'    thresholded spiketimes (spikes where speed is > 5 cm)
%   'SpikeTrain'    original spiketrain

SpikeTrain_thresh = cell(1,length(position));
for sessNum = 1:length(position)
    if ~isempty(position{1,sessNum}) && ~isempty(SpikeTimes{1,sessNum})
        sessionArray = cell(1,length(SpikeTimes{1,sessNum}));
        for unitNum = 1:length(SpikeTrain{1,sessNum})
            ST_now = SpikeTimes{1,sessNum}{1,unitNum};
            SpikeTrain_thresh{1,sessNum} = binSpikes(position{1,sessNum}(:,1), SpikeTimes{1,sessNum});
        end
    end
end

end

