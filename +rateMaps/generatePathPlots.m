function generatePathPlots(pos_cm, trialType, hwLoc)
%GENERATEPATHPLOTS Summary of this function goes here
%   Detailed explanation goes here

for sessNum = 1:length(pos_cm)
    figTit = strcat('\SESS:', sprintf('%.f', sessNum), '\TYPE:', trialType{1,sessNum}, '\hwLoc:', sprintf('%.f', hwLoc{1,sessNum}));
    plot(pos_cm{1,sessNum}(:,2),pos_cm{1,sessNum}(:,3))
    title(figTit)
    pause
    clf
end

end

