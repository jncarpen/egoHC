% create new subplot
% for animal 24116, im starting after the training sessions (at session 7)

hf = figure('units','normalized','outerposition',[0 0 1 1]);
hf = colordef(hf,'white'); %Set color scheme
hf.Color='w'; %Set background color of figure window
% subplot(5,5,1)

[hwLoc, rdLoc] = getWellLoc(labNotes, trialType); % for all sessions

sessNum = 25;
unitNum = 6;
totalUnits = length(SpikeTrain{1,sessNum});
[binCtrs_20, tcVals_20] = goalDirSar(pos{1,sessNum}, hwLoc{1,sessNum}, hd{1,sessNum}, SpikeTimes{1,sessNum}{1,unitNum}, 20);
[binCtrs_30, tcVals_30] = goalDirSar(pos{1,sessNum}, hwLoc{1,sessNum}, hd{1,sessNum}, SpikeTimes{1,sessNum}{1,unitNum}, 30);
[binCtrs_10, tcVals_10] = goalDirSar(pos{1,sessNum}, hwLoc{1,sessNum}, hd{1,sessNum}, SpikeTimes{1,sessNum}{1,unitNum}, 10);
[binCtrs_50, tcVals_50] = goalDirSar(pos{1,sessNum}, hwLoc{1,sessNum}, hd{1,sessNum}, SpikeTimes{1,sessNum}{1,unitNum}, 50);
[binCtrs_05, tcVals_05] = goalDirSar(pos{1,sessNum}, hwLoc{1,sessNum}, hd{1,sessNum}, SpikeTimes{1,sessNum}{1,unitNum}, 5);
[binCtrs_15, tcVals_15] = goalDirSar(pos{1,sessNum}, hwLoc{1,sessNum}, hd{1,sessNum}, SpikeTimes{1,sessNum}{1,unitNum}, 15);


plot(binCtrs_05, tcVals_05, 'Color',[0.3010, 0.7450, 0.9330], 'LineWidth', 1.15, 'DisplayName','5')
hold on
plot(binCtrs_10, tcVals_10, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.15, 'DisplayName','10')
plot(binCtrs_15, tcVals_15, 'Color', 'm', 'LineWidth', 1.15, 'DisplayName','15')
plot(binCtrs_20, tcVals_20, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.15, 'DisplayName','20')
plot(binCtrs_30, tcVals_30, 'Color',[0.4940, 0.1840, 0.5560], 'LineWidth', 1.15, 'DisplayName','30')
plot(binCtrs_50, tcVals_50, 'Color',[0.4660, 0.6740, 0.1880], 'LineWidth', 1.15, 'DisplayName','50')

hold off
% legend
[hleg,att] = legend('show');
% title(hleg,'num_bins')
% hleg.Title.Visible = 'on';

title("goal direction (head direction)")
xlabel("angle (rad)")
ylabel("fr (Hz)")
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
box off