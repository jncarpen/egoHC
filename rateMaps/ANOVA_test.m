% anova visualization

% create data matrix (1xn array, where n = total # of samples)
data = [c{1,1}, c{2,1}];

% cell array with groups
group = {'G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G1','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2','G2'};

% run ANOVA
[p1] = anova1(data, group,'on'); %// Use the 'off' option to prevent the table/box plot from showing up.

% make notched boxplot
boxplot(data, group, 'notch', 'on') 