% from adrien peyrache's data 
% find a good example head direction cell from the AD thalamus

% make position vector
position_HD = position;
pos_hd = zeros(length(trackingtimes), 5);
pos_hd(:,1) = trackingtimes'/1000; % convert to seconds
pos_hd(:,2:5) = position_HD;
sampleRate = mean(diff(trackingtimes/1000));

% Make pathplot_HD
for unitNum =46:length(cellspikes)
    unitNum
    STNow = cellspikes{1,unitNum}/1000;
    if length(STNow) > 20
        figTit = strcat('A12 U', sprintf('%.f', unitNum));
        fig = figure('Position', [100 100 600 300]);
        set(gcf,'color','w');
        sgtitle(figTit)
        
        subplot(1,2,1)
        pathPlot_HD(pos_hd, STNow', rad2deg(headangle)');
        pbaspect([1 1 1])
        box off
        
        subplot(1,2,2)
        spkang = getSpikeAngle(rad2deg(headangle)', STNow');
        tc_HD = analyses.turningCurve(spkang, pos_hd, sampleRate,'smooth', 3, 'binWidth', 10);
        plot(tc_HD(:,1), tc_HD(:,2), 'Color', 'k', 'LineWidth', 1.10)
        pbaspect([1 1 1])
        xlim([0 360])
        xlabel("head angle (deg)")
        ylabel("fr (hz)")
        box off
        
%         pause
%         clf
        
        filename = strcat('D:\egoAnalysis\refUnits\HD\A12_U', sprintf('%.f', unitNum), '.png');
        saveas(fig, filename);
    end
end

