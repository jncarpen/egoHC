function tc = hdTuning(P,HD,ST, nBins)
%HDTUNING: make HD tuning plot
%   INPUTS
%   hd:                     Tx1 double, where T is the number of
%                           timestamps. Each row corresponds to a head
%                           direction value from 1-360 degrees.
%
%   pos:                    [t x y x1 y1]; each column should be of length
%                           T to match the lenght of the hd and SpikeTrain
%                           vectors.
%
%   SpikeTrain:             SpikeTimes that were binned by the number of
%                           position samples- should be a Tx1 double.
%   OUTPUTS
%   TC:                     Head direction tuning curve will be generated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unpack inputs
t = P(:,1);
tpf = mode(diff(t));

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
spk_hd = HD(spkidx);

%% tuning curve
angEdges = linspace(0,360,nBins+1);
[bincounts, edges, binidx] = histcounts(spk_hd,angEdges);

% find standard deviation in each bin
for bin = 1:length(edges)-1
    stdev_spk(bin) = std(inst_fr(find(binidx==bin)), 'omitnan');
end
% occupancy map
for bin = 1:length(angEdges)-1
    idx_now = find(HD>angEdges(bin) & HD<=angEdges(bin+1));
    stdev_bin(bin) = std(inst_fr(idx_now), 'omitnan');
    occ(bin) = sum(P_logical(idx_now)); % occupancy (counts)
end
% bin centers
for i = 1:length(edges)
    if i+1 <= length(edges)
        ctrs(i) = ((edges(i+1)-edges(i))/2)+edges(i);
    end
end
% tuning curve
tc = bincounts./(occ*tpf + eps); 
tc = imgaussfilt(tc, 2, 'Padding', 'circular');

% %% plot
% figure; hold on;
% set(gcf,'color','w');
% errorbar(ctrs, tc, stdev_spk, 'Color', 'k', 'LineWidth', .75);
% plot(ctrs, tc, 'Color', 'k', 'LineWidth', .75);
% ylabel("fr (Hz)"); xlabel("angle (deg)");
% xticks([90 180 270 360]); 
% box off;
% set(gca, 'FontSize', 15, 'FontName', 'Myriad Pro');
end

