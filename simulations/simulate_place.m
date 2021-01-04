function [sim] = simulate_place(param)
%SIMULATE_PLACE: simulate a place cell.
%   INPUTS - 
%   param.position:         [t x y x2 y2]
%   param.ctr_mass:         center of place field, [x,y]
%   param.noise:            proportion of noise to add; can range
%                           from 0 to 1
%   param.width:            width of the place field
%
%   OUTPUTS - 
%   sim.spiketimes:         spiketimes of simulated cells, in seconds
%   sim.position:           position stamps
%
%   Jo Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pull out position information
t = param.position(:,1);
fs = mode(diff(t));
xpos = param.position(:,2);
ypos = param.position(:,3);

% pull out mu (place field center of mass)
xmu = param.ctr_mass(1,1);
ymu = param.ctr_mass(1,2);

% bin the arena and find bin centers
nBins = 100;
xEdges = linspace(0,150,nBins+1);
yEdges = linspace(0,150,nBins+1);
xCenter = (diff(xEdges)./2)+xEdges(1:end-1);
yCenter = (diff(yEdges)./2)+yEdges(1:end-1);
[X,Y] = meshgrid(xCenter,yCenter);

% find bin for each time point
[occ,~,~,binX,binY] = histcounts2(xpos,ypos,xEdges,yEdges);

xvar = nanvar(xpos); % var on x axis
yvar = nanvar(ypos); % var on y axis

xsd = sqrt(xvar); % std deviation on x axis
ysd = sqrt(yvar); % std deviation on y axis

% kendall correlation (excludes NaNs)
% rho_kendall = corr(xpos,ypos,'Type','Kendall','Rows','complete');
rho = corr(xpos,ypos,'Type','Pearson','Rows','complete');
 if (abs(rho) >= 1.0)
        disp("error: rho must lie between -1 and 1");
    return
 end
 
covxy = rho*xsd*ysd; % calculation of the covariance
scaleCov = .75; % scale the covariance matrix by some factor
C = [xvar covxy; covxy yvar].*scaleCov; % the covariance matrix
A = inv(C); % the inverse covariance matrix

% Compute value of Gaussian pdf at each point in the grid
clear zMapped
width = param.width; % values over 1 will give smaller place fields
z = 1/(2*pi*sqrt(det(C))) * exp(-width * (A(1,1)*(X-xmu).^2 + 2*A(1,2)*(X-xmu).*(Y-ymu) + A(2,2)*(Y-ymu).^2));
[zMapped] = map0to1(z); % dependency

clear spiketimes
n = zeros(length(t),1);
n_map = zeros(nBins,nBins);
count = 1;
for frame = 1:length(t)
    binX_now = binX(frame); binY_now = binY(frame);
    if binX_now > 0 && binY_now > 0 
        ri = rand(1);
        pdx_now = zMapped(binX_now, binY_now);
        if ri < pdx_now
            n(frame) = 1;
            n_map(binX_now, binY_now) = n_map(binX_now, binY_now) + n(frame);
            foreground_spikes(count) = t(frame);
            count = count + 1;
        else 
            n(frame) = 0;
        end
    else
        n(frame) = NaN;
    end
end
fr_map = n_map./(occ.*fs);


%% add some noise
throw_away = 1-param.noise;
sz = floor(length(t)-length(t)*throw_away);
background_spikes = datasample(t, sz, 'Replace', false);
% add them together and sort
st_unsorted = [foreground_spikes'; background_spikes];
spiketimes = sort(st_unsorted, 'ascend'); % simulated timestamps


%% visualize
figure; set(gcf,'color','w');
pathPlot_hd(param.position, spiketimes, get_hd(param.position))
title("") % remove title

figure; set(gcf,'color','w');
map = analyses.map(param.position, spiketimes, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
c = colorbar; c.FontName = 'Helvetica'; c.FontSize = 15;
colormap(gca,'jet')
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Helvetica', 'FontSize', 15, 'FontWeight', 'normal');
box off

% subplot(2,2,3)
% surf(X,Y,fr_map,'FaceAlpha',0.5);
% pbaspect([1 1 1])

%% save stuff
sim.spiketimes = spiketimes;
sim.position = param.position;

end 
