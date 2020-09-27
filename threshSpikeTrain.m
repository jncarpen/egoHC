SpikeTrain_thresh = cell(1,length(pos_cm));

for sessNum = 1:length(pos_cm)
    if ~isempty(pos{1,sessNum}) && ~isempty(SpikeTimes_thresh{1,sessNum})
        sessionArray = cell(1,length(SpikeTimes{1,sessNum}));
        for unitNum = 1:length(SpikeTrain{1,sessNum})
            ST_now = SpikeTimes_thresh{1,sessNum}{1,unitNum};
            SpikeTrain_thresh{1,sessNum} = binSpikes(pos_cm{1,sessNum}(:,1), SpikeTimes_thresh{1,sessNum});
        end
    end
end
