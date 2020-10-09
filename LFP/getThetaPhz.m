function [thetaPhz, spkThetaPhz, tc_TP] = getThetaPhz(rawEEG, SpikeTimes)

%THETAPHZ Filter signal within theta range(4-12 Hz).

%   USAGE           thetaFiltSig, phases] = thetaFiltered(lfp)
%
%   INPUTS
%   lfp:            LFP signal (unfiltered) for a SINGLE channel.
%
%   OUTPUTS
%   filtered:       Filtered LFP signal
%
%   Jordan Carpenter
%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%

% add BNT to path (to compute TC)
addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));

% rename input (trying to fix a bug)
lfp = rawEEG;

% parse relevant information
lfp_sig = lfp(:,1); % LFP signal
lfp_t = lfp(:,2); % timestamps
lfp_fr = mean(diff(lfp_t)); % sample rate (diff wasnt working)

% set values (in Hz)
nyquistLim = 125; % for Axona
filterAbv = 4;
filterBlw = 12; 

% filter the signal with a 3rd order bandpass butterworth filter
[b, a]= butter(3,[filterAbv/nyquistLim  filterBlw/nyquistLim],'bandpass');
thetaFiltSig = filtfilt(b, a, lfp_sig);

% grab the phase of the filtered signal at each timestamp
phases = angle(hilbert(thetaFiltSig));

% get theta phase for each time the cell spikes
idx = knnsearch(lfp_t, SpikeTimes);
spkThetaPhz = phases(idx);

% catch errors from knnsearch
if length(SpikeTimes) ~= length(spkThetaPhz)
    disp("Error: SpikeTimes and resulting spkThetaPhz vector are different lengths.")
end

% make thetaPhz matrix [thetaPhz, timestamp]
thetaPhz = [phases, lfp_t];

% make tuning curve
nBins = 20;
edges = linspace(-pi,pi,nBins);
[spkTPmap, mapAxis] = histcounts(spkThetaPhz,edges);
[allPhasesMap] = histcounts(phases,edges);

% compute bin centers
    for i = 1:length(mapAxis)
        if i+1 <= length(mapAxis)
            binCenters(i) = ((mapAxis(i+1)-mapAxis(i))/2)+mapAxis(i);
        end
    end
    
% assign values to tc matrix
tc_TP(:,1) = binCenters;
tc_TP(:,2) = spkTPmap./(allPhasesMap*lfp_fr + eps);
end

