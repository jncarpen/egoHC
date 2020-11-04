function [filtered] = ripFiltered(lfp)

%RIPFILTERED Filter signal between 100-200 Hz.

% Nyquist = 125 (with sampling rate = 250 Hz)
% since the sampling rate for axona is too low (every 0.0040 s), we make a high-pass filter
% to return signals greater than 100 Hz. The idea is to filter within the
% 100-200 Hz range, this can be changed with Neuralynx if there is a
% different sampling frequency. 
% will need to make the butterworth a bandpass filter for signals with
% higher sampling frequency. 
% Add timestamps as an input to this function so that we can calculate the
% sampling frequency and Nyquist limit for each EEG channel that is input.

%   USAGE:
%
%   INPUT:
%   lfp:            LFP signal (unfiltered) for a SINGLE channel.

%   OUTPUT:
%   filtered:       Filtered LFP signal

%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%


nyquistLim = 125; 
filterAbv = 100;

[b, a]= butter(3,[filterAbv/nyquistLim],'high'); % make a butterworth filt
filtered = filtfilt(b, a, lfp);


end

% % for a bandpass filter
% inputsig = eeg_voltage;
% order = 3;
% fcutlow=45;
% fcuthigh=55;
% 
% [b,a]=butter(order,[fcutlow,fcuthigh]/nyquistLim,'bandpass');
% filtsig=filter(b,a,inputsig);  %filtered signal



