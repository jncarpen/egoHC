function [I_rate, I_content] = spatial_info(P, ST)
%SPATIAL_INFO Calculate the *spatial* information rate and content of a
% unit. Adapted from Skaggs et al., 1993
%   INPUTS:
%   'P'                     position vector [t x1 y1 x2 y2]
%   'ST'                    spike times (s)
%   OUTPUTS:
%   'I_rate'                bits per second
%   'I_content'             bits per spike
%   DEPENDENCIES:
%   BNT (Moser Group)       -- @todo remove dependencies
%   J. Carpenter, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = P;
spiketimes = ST;

% calculate ratemap
nbins = 30; % 5 cm bins (900 total bins)
sigma = 2; % smoothing factor
map = analyses.map(position, spiketimes,'binWidth', 1.50/nbins); % bug in 2d smoothing fn
spikemap = imgaussfilt(map.Nspikes, sigma); % number of spikes per spatial bin (smoothed)

% get rid of NaNs in the occupancy map
occ_raw = map.timeRaw;
occ_raw(isnan(occ_raw))=0;
occupancy = imgaussfilt(occ_raw, sigma);

% calculate probability density for rat being at each loc
pdx = occupancy ./ sum(occupancy, 'all', 'omitnan');

% calculate lambda_x (mean fr when rat is at loc x)
lambda_x = eps + (spikemap./(occupancy + eps));

% calculate [overall] mean firing rate (lambda) of the unit
lambda = sum(spikemap, 'all', 'omitnan')./sum(occupancy, 'all', 'omitnan');

% calculate information rate of the cell (bits/s)
I_rate = sum(lambda_x .* log2(lambda_x ./ lambda) .* pdx, 'all', 'omitnan');

% calculate information content (bits/spk)
I_content = I_rate./lambda;

end

