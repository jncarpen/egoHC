function vizEEG(rawEEG, sessNum)

%VIZEEG Visualize all EEG channels at an offset, 4 seconds at a time.
%   Click through plots of all EEG channels at once
%   This may help in finding a channel to use for ripple detection by
%   visualizing which channel has the highest amplitude ripples.

%   INPUT:
%   rawEEG:     Cell array with each channel [rawSignal timeStamp]
%   sessNum:    Which sessNum do you want to look at?

%   OUTPUT:
%   plots:      EEG signals offset every 1000 ms (EEG channels are plotted
%               with the lowest number (i.e. EEG channel in position 1) at
%               the bottom of the window while EEG channel with the highest
%               number is plotted at the top.

%   At the moment, this script works best for the way I've read in Axona
%   data, for other recording softwares, this may need editing.

%   ADD:    raster plot with spike times underneath the traces ** (to look for
%           population bursts which may help you look for ripples visually.

%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%


numChan = length(rawEEG{1,sessNum}); % grab # of channels

    for ts = 1:1000:length(rawEEG{1,sessNum}{1,1})
        for i=1:numChan
            plot(rawEEG{1,sessNum}{1,1}(ts:ts+1000,2),100*i+rawEEG{1,sessNum}{1,i}(ts:ts+1000,1),'k')
            title("vizEEG")
            xlabel("time (s)")
            ylabel("amplitude")
            hold on
        end
        pause
        clf
    end
end


% find sample rate: 1000./diff(rawEEG{1,13}{1,1}(1:10,2))
