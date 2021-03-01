function [tc, ctrs] = plot_egoBearing(P, ST, rp, doPlot)
%PLOT_EGOBEARING Thresholded!
% October 22, 2020
% This function can probably replace the other one (egoBearing.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unpack inputs
t = P(:,1); tpf = mode(diff(t));
X = P(:,2); Y = P(:,3);
rlX = rp(1); rlY = rp(2);

% head direction
[HD] = get_hd(P); % in deg (0 is 'up')

% allo & egocentric bearing
allo = mod(atan2d(rlY-Y, rlX-X), 360);
allo = mod(allo + 90, 360); % shift so that 0 deg is in front of animal
ego = mod(allo - HD, 360);

% spiketrain
startT = t(1); stopT = t(end);
ST = ST(ST < stopT & ST > startT);
time_edges = linspace(startT,stopT,numel(t)+1);
spiketrain = histcounts(ST, time_edges);

% instantaneous fr at each timepoint
inst_fr = spiketrain./tpf;

%% apply speed threshold 
% @todo: make this an option
threshold = 5; % cm/s
[P_logical, ST_thresh] = speed_threshold(P, ST, threshold);

% find time indices when cell spikes
spkidx = knnsearch(t, ST_thresh);
spk_ego = ego(spkidx);

%% compute tuning curve
nBins = 40; % 9 degree bins
angEdges = linspace(0,360,nBins);
[bincounts, edges, bin] = histcounts(spk_ego,angEdges);

% find standard deviation in each bin
for binNum = 1:length(edges)-1
    stdev_spk(binNum) = std(inst_fr(find(bin==binNum)), 'omitnan');
end
% calculate occupancy map
for bin = 1:length(angEdges)-1
    idx_now = find(ego>angEdges(bin) & ego<=angEdges(bin+1));
    stdev_bin(bin) = std(inst_fr(idx_now), 'omitnan');
    ego_occ(bin) = sum(P_logical(idx_now)); % occupancy (counts)
end
for i = 1:length(edges)
    if i+1 <= length(edges)
        ctrs(i) = ((edges(i+1)-edges(i))/2)+edges(i);
    end
end
% tuning curve
tc = bincounts./(ego_occ*tpf + eps); 
tc = imgaussfilt(tc, 2, 'Padding', 'circular');

%% plot

switch doPlot
    case "True"
        hold on;
        plot(ctrs, tc, 'Color', [0 .4 .9], 'LineWidth', 1.10)
%         errorbar(tc,stdev_spk);
        ylabel("fr (Hz)"); xlabel("angle (deg)");
        xticks([90 180 270 360])
        box off
    case "False"
        % do nothing
end


end

