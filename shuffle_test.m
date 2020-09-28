sessNum =27; unitNum = 1;
STNow = SpikeTimes_thresh{1,sessNum}{1,unitNum};

totalIter = 30; shiftVal = 150;
for iter = 1%:totalIter
    r = STNow + 30;
    shiftVal = shiftVal + (shiftVal/totalIter);
end