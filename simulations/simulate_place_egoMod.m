function [sim] = simulate_place_egoMod(param)
%SIMULATE_PLACE_EGO: simulate a place cell.
%   INPUTS - 
%   param.position:         [t x y x2 y2]
%   param.ctr_mass:         center of place field, [x,y]
%   param.noise:            proportion of noise to add; can range
%                           from 0 to 1
%   param.width:            width of the place field
%   param.theta:            preferred head-direction
%
%   OUTPUTS - 
%   sim.spiketimes:         spiketimes of simulated cells, in seconds
%   sim.position:           position stamps
%
%   Jo Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% turn off annoying warnings
warning('off','all')

% pull out position information
t = param.position(:,1);
fs = mode(diff(t));
xpos = param.position(:,2);
ypos = param.position(:,3);

% pull out mu (place field center of mass)
xmu = param.ctr_mass(1,1);
ymu = param.ctr_mass(1,2);

% bin the arena and find bin centers
nBins = 50;
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

% kendall correlation (excludes NaNs) @joCarp (kinda slow)
rho = corr(xpos,ypos,'Type','Pearson','Rows','complete');
 if (abs(rho) >= 1.0)
        disp("error: rho must lie between -1 and 1");
    return
 end
 
covxy = rho*xsd*ysd; % calculation of the covariance
scaleCov = 1; % scale the covariance matrix by some factor
C = [xvar covxy; covxy yvar]./scaleCov; % the covariance matrix
A = inv(C); % the inverse covariance matrix

% Compute value of Gaussian pdf at each point in the grid
width = param.width; % values over 1 will give smaller place fields
z = 1/(2*pi*sqrt(det(C))) * exp(-width * (A(1,1)*(X-xmu).^2 + 2*A(1,2)*(X-xmu).*(Y-ymu) + A(2,2)*(Y-ymu).^2));
[zMapped] = map0to1(z); 

n = zeros(length(t),1);
n_map = zeros(nBins,nBins);
count = 1;
clear spiketimes
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

% calculate egocentric bearing at each spiketime
head_direction = get_hd(param.position); % shift the values 
idx = knnsearch(t,foreground_spikes');
spkhd = head_direction(idx);


%% find the probability of each hd-angle from the von-Mises distribution
theta_pref = deg2rad(param.theta); % vm takes radian inputs
kappa = 5;
[spkhd_pdf, ~] = circ_vmpdf(deg2rad(spkhd), theta_pref, kappa); % circstat toolbox
[pdf_mapped] = map0to1(spkhd_pdf); 
fspikes_transposed = foreground_spikes';

[B,Bsort]=sort(deg2rad(spkhd), 'ascend');
C = spkhd_pdf(Bsort);
figure; hold on;
plot(rad2deg(B),C, '.', 'Color', [.5 .5 .5])
xticks([90 180 270 360])
xline(param.theta, ':r','LineWidth', 2);
hold off;

clear foreground_spikes_saved
count = 1;
for h = 1:length(spkhd)
    ri = rand(1);
    if ri < pdf_mapped(h)
        fspikes(count) = fspikes_transposed(h);
        count = count + 1;
    end
    
end


%% add some noise
throw_away = 1-param.noise;
sz = floor(length(t)-length(t)*throw_away);
background_spikes = datasample(t, sz, 'Replace', false);
% add them together and sort
st_unsorted = [fspikes'; background_spikes];
spiketimes = sort(st_unsorted, 'ascend'); % simulated timestamps


%% make a [speed-filtered] spike train

% speed threshold spikes (keep >5cm/s)
x_smooth=medfilt1(xpos); y_smooth=medfilt1(ypos);
speed = zeros(length(t), 1);
for i = 2:numel(x_smooth)-1
    speed(i) = sqrt((x_smooth(i+1) - x_smooth(i-1))^2 + (y_smooth(i+1) - y_smooth(i-1))^2) / (t(i+1) - t(i-1));
end

stopTime = nanmax(t); startTime = nanmin(t);
SpikeTimes_sim = spiketimes(spiketimes < stopTime & spiketimes > startTime); % remove times outside of recording
edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate
binnedSpikes = histcounts(SpikeTimes_sim, edgesT); % bin those spikes baby
speed_idx = find(speed<5); % find indices when animal was moving slow
binnedSpikes(speed_idx)=0; % get rid of spikes when animal was moving slow
binnedSpikes = imgaussfilt(binnedSpikes, 2, 'Padding', 'replicate'); % smooth ST
SpikeTrain_sim = binnedSpikes;

% seperate spiketimes
% ST_pref = 
% ST_nonpref = 


%% visualize (with colormap)

figure; set(gcf,'color','w');
figTit = strcat('preferred theta =', {' '},  sprintf('%.f', param.theta), {' '}, 'degrees');
% sgtitle(figTit{1,1}, 'FontName', 'Calibri light', 'FontSize', 25, 'FontWeight', 'normal')

subplot(1,2,1)
hold on;
pathPlot_hd(param.position, spiketimes, get_hd(param.position))
title("")
hold off;

subplot(1,2,2)
map = analyses.map(param.position, spiketimes, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
colorbar
colormap(gca,'jet')
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
box off

%% visualize (with arrow spike plot)
% figure; set(gcf,'color','w');
% figTit = strcat('preferred theta =', {' '},  sprintf('%.f', param.theta), {' '}, 'degrees');
% % sgtitle(figTit{1,1}, 'FontName', 'Calibri light', 'FontSize', 25, 'FontWeight', 'normal')
% 
% subplot(1,2,1)
% hold on;
% pathPlot_quiver(param.position, spiketimes, get_hd(param.position))
% hold off;
% 
% subplot(1,2,2)
% map = analyses.map(param.position, spiketimes, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
% peakRate = nanmax(nanmax(map.z));
% rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
% plot.colorMap(map.z)
% pbaspect([1 1 1])
% colorbar
% colormap(gca,'jet')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
% box off


%% save stuff
sim.spiketimes_unfiltered = spiketimes;
sim.spiketimes = SpikeTimes_sim;
sim.spiketrain = SpikeTrain_sim;
sim.position = param.position;
% sim.hd = hd_sim;

% turn back on warnings
warning('on','all')

end 


%% SCRATCH
% % plot the pdf
% figure;
% obj1 = CircHist(spkhd, 100);
% figure;
% x = 1:0.02:360;
% mu = param.theta; sigma = 45;
% y = normpdf(x,mu,sigma);
% plot(x,y)