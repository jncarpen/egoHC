function hdTuning(hd, pos, SpikeTrain)
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

% make sure that everything is in SECONDS ** for appropriate Hz
% calculations :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert degrees to radians 
if max(hd) > 2*pi
    hd = hd*pi/180;
end

trackingtimes = pos(:,1); % pull trackingtimes out of pos matrix

% find video sampling frequency (in S)
sampleRate = mode(diff(trackingtimes));

%angular bins
da = pi/30; %6 degrees
angBins = [da/2:da:2*pi-da/ 2];

%Occupancy
[histAng,~,~] = histcounts(hd,angBins+1);

% histogram of the number of times the cell fired in each bin of
% head-direction
[spkPerAng,~,~] = histcounts(SpikeTrain,angBins);

% now compute the tuning curve:
HDT = spkPerAng./histAng*sampleRate;

angBins = angBins(2:end); % cut out the first bin

figure(1),clf
%set(gcf,'Position',[62   319   783   281])
subplot(1,3,1)
    plot(angBins,spkPerAng)
% polar(angBins,spkPerAng)
    title('Number of spikes')
subplot(1,3,2)
    plot(angBins,histAng)
% polar(angBins,histAng)
    title('Occupancy')
subplot(1,3,3)
    plot(angBins,HDT)
% polar(angBins,hdTuning)
    title('Tuning Curve (Hz)')
return

end

