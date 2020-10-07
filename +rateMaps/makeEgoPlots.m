% need to specify reference locations

for sessNum = 1:length(SpikeTimes)
    ST = SpikeTimes{1,sessNum};
     if ~isempty(ST) % skip empty trials
         for unitNum = 1:length(ST)
            fig = egoPlots(pos,SpikeTimes, refLoc, hd, sessInfo, sessNum, unitNum);
             
            % name figure
            UID = UniqueID{1,sessNum}{1,unit}; 
            dateStr = char(extractBetween(sessInfo{1,sessNum}.trial_date{1,1},",","20"));
            dateStr = dateStr(find(~isspace(dateStr)));
            timeStr = sessInfo{1,sessNum}.trial_time{1,1}(1:5);
            figTit = strcat('UID:', sprintf('%.f', UID), '\SESS:', sprintf('%.f', sessNum), '\DAT:', dateStr, '\TIME:', timeStr, '\TYPE:', trialType{1,sessNum});
            fileBody = strcat('UID', sprintf('%.f', UID), '_SESS', sprintf('%.f', sessNum), '_TYPE_', trialType{1,sessNum});
            fig.Name = figTit; % set figure name

            % save figures
            filename = strcat('D:\egoPlots\', fileBody, '.png');
            saveas(fig, filename);

            clf
            close all
             
         end
     end
end
         