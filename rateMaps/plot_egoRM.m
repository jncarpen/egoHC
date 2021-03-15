function out = plot_egoRM(P, ST, RP)
%EGORM 
%   INPUTS-
%   P:      position vector
%   ST:     spiketimes
%   RP:     reference point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
% if RP is output from the cosine model
% convert to cm
% nbinstheta = 10;
% RP = RP.*(150/nbinstheta);

ego_RM = egoRateMap(P, ST, RP);
[Xi,Yi,Zi,Ci] = polarplot3d(ego_RM, 'PlotType','off');

fig = surf (Xi,Yi,Zi,Ci,'LineStyle','none'); hold on;
set(gcf,'color','w');
set(gca, 'visible', 'off', 'box', 'off');
view(2);
camroll(90);
pbaspect([1 1 1]);
grid off
shading interp

% find peak firing rate
maxfr = nanmax(Ci(:));
[row, col] = find(ismember(Ci, maxfr));
locxy = [Xi(row,col), Yi(row,col), Zi(row,col)];
% scatter3(locxy(1), locxy(2), locxy(3), [30], 'k', 'filled');

% package output
out.fig = fig;
out.maxfr = maxfr;
out.loc = locxy;

end

