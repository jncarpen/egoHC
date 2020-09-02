function alignPosition(sessInfo, pos)
%ALIGNPOSITION Summary of this function goes here
%   Detailed explanation goes here
%
%    J. Carpenter 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse sessInfo struct
window_min_x = []; window_max_x = [];
window_min_y = []; window_max_y = [];
min_x = []; max_x = [];
min_y = []; max_y = [];

% pool max/min information across sessions
for sess = 1:length(sessInfo)
    window_min_x = [window_min_x, str2double(sessInfo{1,sess}.window_min_x{:})];
    window_max_x = [window_max_x, str2double(sessInfo{1,sess}.window_max_x{:})];
    window_min_y = [window_min_y, str2double(sessInfo{1,sess}.window_min_y{:})];
    window_max_y = [window_max_y, str2double(sessInfo{1,sess}.window_max_y{:})];

    min_x = [min_x, str2double(sessInfo{1,sess}.min_x{:})];
    max_x = [max_x, str2double(sessInfo{1,sess}.max_x{:})];
    min_y = [min_y, str2double(sessInfo{1,sess}.min_y{:})];
    max_y = [max_y, str2double(sessInfo{1,sess}.max_y{:})];
end

% define mean window positions
mean_win_min_x = floor(mean(window_min_x));
mean_win_max_x = floor(mean(window_max_x));
mean_win_min_y = floor(mean(window_min_y));
mean_win_max_y = floor(mean(window_max_y));

for col = 1:length(window_min_x)
    diff_win_min_x = window_min_x - mean_win_min_x;
    diff_win_max_x = window_max_x - mean_win_max_x;
    diff_win_min_y = window_min_y - mean_win_min_y;
    diff_win_max_y = window_max_y - mean_win_max_y;
end

diff_x = (diff_win_max_x-diff_win_min_x)';
diff_y = (diff_win_max_y-diff_win_min_y)';

for sess = 1:length(pos)
    shft_pos{1,sess}(:,2) =  pos{1,sess}(:,2) + diff_x(sess); % L1A(x)
    shft_pos{1,sess}(:,3) =  pos{1,sess}(:,3) + diff_x(sess); % L1B(y)
end

% plot mean window
figure
i = 25;
xline(mean_win_min_x, 'r', 'LineWidth', 1.25);
hold on
plot(shft_pos{1,i}(:,2), shft_pos{1,i}(:,3))
plot(pos{1,i}(:,2), pos{1,i}(:,3))

xline(mean_win_max_x, 'r', 'LineWidth', 1.25);
yline(mean_win_min_y, 'r', 'LineWidth', 1.25);
yline(mean_win_max_y, 'r', 'LineWidth', 1.25);


% for i = 1:length(pos)
%     plot(shft_pos{1,i}(:,2), shft_pos{1,i}(:,3))
%     legend
%     pause
% end

% legend
title("Mean window")
% xlim([mean_win_min_x-10, mean_win_max_x+10])
% ylim([mean_win_min_y-10, mean_win_max_y+10])
hold off


end

