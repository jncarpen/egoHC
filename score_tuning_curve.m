function [mean_stats, sum_stats] = score_tuning_curve(stats_array)
%SCORE_TUNING_CURVE Summary of this function goes here
%   INPUTS
%   'stats_array'      cell array of tuning curve statistics structs (BNT)
%   'measure'          which measure to compute
%   OUTPUT
%   'mean_stat'         mean of statistic across all 2D spatial bins
%   Jordan Carpenter, October 5, 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate length of input
[numRow, numCol] = size(stats_array); 

% MVL
k = [];
for row = 1:numRow
    for col = 1:numCol
        k = [k, stats_array{row,col}.r];
    end
end
mean_stats.MVL = mean(k, 'omitnan');
sum_stats.MVL = sum(k, 'omitnan');

% mean direction
k = [];
for row = 1:numRow
    for col = 1:numCol
        k = [k, stats_array{row,col}.mean];
    end
end
mean_stats.mean_direction = mean(k, 'omitnan');
sum_stats.mean_direction = sum(k, 'omitnan');

% standard deviation
k = [];
for row = 1:numRow
    for col = 1:numCol
        k = [k, stats_array{row,col}.std];
    end
end
mean_stats.std = mean(k, 'omitnan');
sum_stats.std = sum(k, 'omitnan');

% peak firing rate
k = [];
for row = 1:numRow
    for col = 1:numCol
        k = [k, stats_array{row,col}.peakRate];
    end
end
mean_stats.peak_rate = mean(k, 'omitnan');
sum_stats.peak_rate = sum(k, 'omitnan');

% mean firing rate
k = [];
for row = 1:numRow
    for col = 1:numCol
        k = [k, stats_array{row,col}.meanRate];
    end
end
mean_stats.mean_rate = mean(k, 'omitnan');
sum_stats.mean_rate = sum(k, 'omitnan');

end

