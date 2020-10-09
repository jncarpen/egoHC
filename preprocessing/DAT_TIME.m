date = cell(1,length(sessInfo));
time = cell(1,length(sessInfo));
sessNum = cell(1,length(sessInfo));

for sess = 1:length(sessInfo)
    sessNum{1,sess} = sess;
    date{1,sess} = sessInfo{1,sess}.trial_date{1,1};
    time{1,sess} = sessInfo{1,sess}.trial_time{1,1};
end

SDT = table(sessNum', date', time', {'sessNum', 'Date', 'Time'});