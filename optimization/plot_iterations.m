function plot_iterations(model)
% plot the output of the model at each iteration

% grab struct
g = model.saved;
b = model.bestParams;

% take the mean of the errors
err = mean(g.FV,2);

% make a title
% plot_title = strcat('xref=', {' '}, sprintf('%.f', b.xref(end)), 'yref=', {' '}, sprintf('%.f', b.yref(end)));
plot_title = strcat('Current Function Value = ', {' '}, sprintf('%.f', err(end)));

hold on;
set(gca, 'FontName', 'Calibri Light', 'FontSize', 20);
title(plot_title{1,1}, 'FontName', 'Calibri light', 'FontSize', 25, 'FontWeight', 'normal');
xlabel('iteration','FontName', 'Calibri light', 'FontSize', 20, 'FontWeight', 'normal');
ylabel('function value (error)','FontName', 'Calibri light', 'FontSize', 20, 'FontWeight', 'normal');
s = scatter(g.IterCnt, err, [8], 'd','MarkerFaceColor','k','MarkerEdgeColor','k');
pbaspect([1 1 1])
alpha(s,.25)
box off

end

