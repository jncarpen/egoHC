function hdTuning(head_angle, trackingtimes, cellSpikeTrain)
%HDTUNING: make HD tuning plot

% make sure that everything is in SECONDS ** for appropriate Hz
% calculations :)

% find video sampling frequency (there may be a better way to do this)
sampleRate = mode(diff(trackingtimes));

%angular bins
da = pi/30; %6 degrees
angBins = [da/2:da:2*pi-da/2];

%Occupancy
[histAng,~,~] = histcounts(head_angle,angBins);

% histogram of the number of times the cell fired in each bin of
% head-direction
[spkPerAng,~,~] = histcounts(cellSpikeTrain,angBins);

% now compute the tuning curve:
hdTuning = spkPerAng./histAng * sampleRate;

angBins = angBins(2:end);

figure(1),clf
set(gcf,'Position',[62   319   783   281])
subplot(1,3,1)
    plot(angBins,spkPerAng)
% polar(angBins,spkPerAng)
    title('Number of spikes')
subplot(1,3,2)
    plot(angBins,histAng)
% polar(angBins,histAng)
    title('Occupancy')
subplot(1,3,3)
    plot(angBins,hdTuning)
% polar(angBins,hdTuning)
    title('Tuning Curve (Hz)')
return

end

