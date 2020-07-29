function [rip, yrx, opt]=detectRipples(y, t, fs, opt)
% function to detect ripple events in lfp trace within specified time
% windows of a recording.
% inputs:
%   - fn0: string argument with filename incl path
%       eg : fn0='E:\data\xyz.ncs'
%   - opt: structure array with options:
%             opt.filetype: filetype of raw data file:
%                   can be: 'ncs'  - cheetah (default)
%                           'mat'  - ML files, Keith
%                           'abf2' - abf2 file Randall lab
%                  eg > opt.filetype='ncs';

%             opt.verbose: activate verbose mode, plot every ripple detection
%                  eg > opt.verbose=1;
%
%             opt.minlen: min length of a ripple in seconds
%                  eg > opt.minlen=0.015;
%
%             opt.maxlen: maximum length of a ripple in seconds
%                  eg > opt.maxlen=0.300;
%
%             opt.mingap: minimum gap for ripples to be considered separate events
%                  eg > opt.mingap=0.010;
%
%             opt.minamp: minimum amplitude for a ripple event in uV
%                  eg > opt.minamp=50;
%
%             opt.detectlim: detection limit in signals STD
%                  eg > opt.detectlim=4;

%             opt.startendlim: local limit for for start end times in signal STD
%                   eg >  opt.startendlim=2;
%
%             opt.load=0; if the script has run once before, it will save
%             the filtered downsampled trace in the source directory, ste
%             thsi option to 1 for loading teh previously filtered and
%             downsampled trace.
%                   eg > opt.load=0;
%
%             opt.reFs: target sampling Frequency for resampling
%                   eg > opt.reFs=1000:
%
%             opt.plt: set to 1 for displaying summary plot of results
%                   eg > opt.plt=1;
%
%             opt.ts: n*2 matrix of time ranges in SECONDS during which spindles are
%             to be detected. N.B. this field is mandatory!
%                   eg > opt.ts = [2340 2400; 2460 3004; 3418 3697];
%
%             opt.save: save filtered downsampled data (1 to save, 0 not to
%             save)
%
%             opt.saveres: save spindle results mat file (1 to save, 0 not
%             to save)

%
% Ullrich Bartsch, University of Bristol, 2010
% Subsequent tweaks by Rich
%
%--------------------------------------------------------------------------
%some settings:

if nargin < 2, opt = struct(); end

pltwin=0.5;
% band=opt.filterband;
%check filtered ripples are smaller than
noise_lim=50;


%check options
optstrings={'verbose','minlen','mingap',...
    'maxlen', 'detectlim','startendlim', ...
    'reFs','plt','minamp', 'filterband', ...
    'tsmooth', 'resfile'};
standards={0,0.015,0.100,0.500,3.5,1.5,1000,1,50,[125, 250], 0.002, 'ripples'};


for i=1:length(optstrings)
    if ~(isfield(opt ,optstrings{i}))
        opt = setfield(opt, optstrings{i} , standards{i});
        disp([ '   > field opt.' optstrings{i} ' not set! Using standard value: ' num2str(standards{i})])
    end
end

%--------------------------------------------------------------------------
% prepare data
[yr,tr,fsr]=resample_eeg(y, t, fs, opt.reFs);

%rectify (using modified version below) -RG
[yrx]=rectify_eeg(yr,tr,opt.filterband,opt.ts);

% RG: replace spline envelope fit with smoothed hilbert amplitude
amp = abs(hilbert(yrx.fz));
nsm = ceil(fsr*opt.tsmooth);
mx2 = rg.signal.gsmooth(amp, [], nsm);
mxU = (mx2 + yrx.fzmean) * yrx.fzstd; % convert back to uV

opt.detectlim0 = opt.detectlim*yrx.fzstd + yrx.fzmean;
opt.startendlim0 = opt.detectlim*yrx.fzstd + yrx.fzmean;

drawnow;

%--------------------------------------------------------------------------
% ripple detection

% detect threshold crossings:
start_ind=find(mx2(1:end-1)<opt.detectlim & mx2(2:end)>opt.detectlim);
stop_ind=find(mx2(1:end-1)>opt.detectlim & mx2(2:end)<opt.detectlim);
if stop_ind(1)< start_ind(1); stop_ind(1)=[]; end

%%%%%%%%%%%%%%%%%%
% ADDED RG
if length(start_ind) > length(stop_ind)
    start_ind(end) = [];
end

% Only use events that are both inside time range
v1 = rg.helpers.isinrange(tr(start_ind), opt.ts);
v2 = rg.helpers.isinrange(tr(stop_ind), opt.ts);
start_ind = start_ind(v1&v2);
stop_ind  = stop_ind(v1&v2);

% /ADDED RG
%%%%%%%%%%%%%%%%%%

if opt.verbose
    %     close (mh)
    fh2 = figure;
    
    %     plot(cts,cscx.rz,'k')
    %     hold on
    %     plot(cts,cscx.z+30,'color',[.9 .9 .9]);
    %     plot(cts,cscx.fz+20,'color',[.7 .7 .7]);
    %     plot(cts,ones(size(cts))*23.5,'r')
    %     plot(cts,mx2,'r.');
    %     plot(cts(start_ind),mx2(start_ind),'go');
    %     plot(cts(stop_ind),mx2(stop_ind),'bo');
    %     hold off
    disp('##########################################################################')
    disp('                   Ripple detection: Verbose mode!'  )
    disp('##########################################################################')
    disp(' ')
    disp('                    Hit any key to step through!')
    disp(' ')
end

% for every detection apply minimum length and minimum gap restrictions
j=1;i=1;
max_n=length(start_ind);
tmp=zeros(size(length(start_ind),2));

while i<length(start_ind)-1
    rip_len=tr(stop_ind(i))-tr(start_ind(i));
    rip_gap=tr(start_ind(i+1))-tr(stop_ind(i));
    
    
    if opt.verbose
        clf
        
        twin = [tr(start_ind(i))-pltwin-1 tr(stop_ind(i))+pltwin+1];
        idx = m_lookup(twin,tr');
        idx(idx == 0) = 1;
        idx = idx(1):idx(2);
        
        plot(tr(idx),yrx.rz(idx),'k')
        hold on
        plot(tr(idx),yrx.z(idx)+30,'color',[.9 .9 .9]);
        plot(tr(idx),yrx.fz(idx)+20,'color',[.7 .7 .7]);
        plot(tr(idx),ones(size(tr(idx)))*23.5,'r')
        plot(tr(idx),mx2(idx),'r.');
        plot(tr(start_ind(i)),mx2(start_ind(i)),'go');
        plot(tr(stop_ind(i)),mx2(stop_ind(i)),'bo');
        
        disp (['crossing ' num2str([ i  start_ind(i) stop_ind(i) rip_len rip_gap    ] ) ]);
        plot(tr(start_ind(i):stop_ind(i)),yrx.rz(start_ind(i):stop_ind(i)),'r');
        hold off
        axis ([tr(start_ind(i))-pltwin,tr(stop_ind(i))+pltwin,-2 40])
        pause
    end
    
    ii=0;
    detect=0;
    
    if ((tr(stop_ind(i))-tr(start_ind(i))) >=opt.minlen) && ((tr(stop_ind(i))-tr(start_ind(i))) <=opt.maxlen)% event is longer than min length?
        tmp(j,1)=start_ind(i);
        tmp(j,2)=stop_ind(i);
        if opt.verbose; disp([  '   > minlength! '] ); end
        
        
        while  (tr(start_ind(i+ii+1)) - tr(stop_ind(i+ii)) <= opt.mingap) && (i+ii+1 <= max_n)  % other events close?
            
            if opt.verbose
                disp( [ '   < mingap! ' ] )
                disp([ '         '  num2str([ i+ii+1 start_ind(i+ii+1) stop_ind(i+ii+1)...
                    tr(stop_ind(i+ii+1))-tr(start_ind(i+ii+1)) ...
                    tr(stop_ind(i+ii+1))-tr(start_ind(i+ii+1))  ]) ] )
            end
            
            tmp(j,2)=stop_ind(i+ii+1);
            ii=ii+1;
            if (i+ii+1)>length(start_ind)
                break
            end
        end
        
        %compute simple max min amplitude
        [mx0,mxk0]=max( yrx.fz(tmp(j,1):tmp(j,2)) );
        [mn0,mnk0]=min( yrx.fz(tmp(j,1):tmp(j,2)) );
        amp0=sqrt(( mx0 - mn0 )^2 );
        if  (amp0*yrx.fzstd>opt.minamp) && (amp0<50)
            
            i=i+ii+1;
            j=j+1;
            detect=1;
        else
            tmp(j,:)=[];
            i=i+1;
            
        end
        
        
    elseif (tr(start_ind(i+ii+1)) - tr(stop_ind(i+ii)) <= opt.mingap) && (i+ii+1 <= max_n) % event is too short but next event is close by
        if opt.verbose ; disp([ '   < mingap! '  ]); end
        tmp(j,1)=start_ind(i);
        tmp(j,2)=stop_ind(i);
        
        while  (tr(start_ind(i+ii+1)) - tr(stop_ind(i+ii)) < opt.mingap) && (i+ii+1 < max_n-2)  % more than one event close by?
            if opt.verbose
                disp([ '   < mingap! '  ]);
                disp([ '         '  num2str([ i+ii+1 start_ind(i+ii+1) stop_ind(i+ii+1)...
                    tr(stop_ind(i+ii+1))-tr(start_ind(i+ii+1))...
                    tr(start_ind(i+ii+1))-tr(stop_ind(i+ii))]) ] );
            end
            
            tmp(j,2)=stop_ind(i+ii+1);
            ii=ii+1;
        end
        
        
        [mx0,mxk0]=max( yrx.fz(tmp(j,1):tmp(j,2)) );
        [mn0,mnk0]=min( yrx.fz(tmp(j,1):tmp(j,2)) );
        amp0=sqrt(( mx0 - mn0 )^2 );
        
        if    ((tr(tmp(j,2)) - tr(tmp(j,1)))  >= opt.minlen)...
                && ((tr(tmp(j,2)) - tr(tmp(j,1))) <=opt.maxlen)...
                && (amp0*yrx.fzstd>opt.minamp)...
                &&  (amp0<noise_lim) % is it altogether long enough. etc?
            j=j+1;
            i=i+ii+1;
            detect=1;
            
        else
            
            tmp(j,:)=[];
            i=i+1;
        end
        %compute simple max min amplitude
        
        
        
        
    else
        i=i+1;
    end
    
    
    
    
    
    if detect
        jj=j-1;
        %local threshold
        k1=[];
        if tmp(jj,1)>100 %min 100 samples from start
            k1=find( mx2(tmp(jj,1)-100:tmp(jj,1)-1)<opt.startendlim &  mx2(tmp(jj,1)-99:tmp(jj,1))>opt.startendlim, 1,'last' );
        end
        if ~isempty(k1)
            tmp2(jj,1)=tmp(jj,1)-(100-max(k1));
        else
            tmp2(jj,1)=tmp(jj,1);
        end
        
        k2=find( mx2(tmp(jj,2):tmp(jj,2)+99)>opt.startendlim &  mx2(tmp(jj,2)+1:tmp(jj,2)+100)<opt.startendlim ,1,'first');
        if ~isempty(k2)
            tmp2(jj,2)=min(k2)-1+tmp(jj,2);
        else
            tmp2(jj,2)=tmp(jj,2);
        end
        
        % closest local minimum?
        
        % check if new start/stop times do not include other detections?
        
        % all local maxima & minima
        [mx0,mxk0]=findpeaks(  yrx.f(tmp2(jj,1):tmp2(jj,2)) );
        [mn0,mnk0]=findpeaks( -yrx.f(tmp2(jj,1):tmp2(jj,2)) );
        [amx0,amxk0]=max( mx2(tmp2(jj,1):tmp2(jj,2)) );
        
        rip(jj).indstart=tmp2(jj,1);
        rip(jj).indstop=tmp2(jj,2);
        
        rip(jj).indallmax=tmp2(jj,1)+mxk0-1;
        rip(jj).indallmin=tmp2(jj,1)+mnk0-1;
        rip(jj).indmaxamp=tmp2(jj,1)+amxk0-1;
        
        [mmx,mmxk]=max(mx0);
        [mmn,mmnk]=max(mn0);
        
        % absoulute local maximum & minimum
        rip(jj).indabsmax=tmp2(jj,1)+mxk0(mmxk)-1;
        rip(jj).indabsmin=tmp2(jj,1)+mnk0(mmnk)-1;
        
        rip(jj).tsstart=tr(rip(jj).indstart);
        rip(jj).tsstop=tr(rip(jj).indstop);
        
        rip(jj).tsabsmax=tr(rip(jj).indabsmax);
        rip(jj).tsabsmin=tr(rip(jj).indabsmin);
        
        rip(jj).tsallmax=tr(rip(jj).indallmax);
        rip(jj).tsallmin=tr(rip(jj).indallmin);
        rip(jj).tsmaxamp=tr(rip(jj).indmaxamp);
        
        rip(jj).amplitude=mxU(rip(jj).indmaxamp);
        %         rip(jj).amplitude=abs (cscx.f(rip(jj).indabsmax))+ abs( cscx.f(rip(jj).indabsmin));
        rip(jj).length=tr(rip(jj).indstop)-tr(rip(jj).indstart);
        rip(jj).freq=length(rip(jj).indallmax)/rip(jj).length;
        
        if opt.verbose==1
            
            figure(fh2)
            hold on
            
            
            axis ([tr(tmp2(jj,1))-pltwin,tr(tmp2(jj,2))+pltwin,-2 40])
            
            plot(tr(rip(jj).indstart:rip(jj).indstop), yrx.rz(rip(jj).indstart : rip(jj).indstop) ,'b');
            plot(tr(rip(jj).indstart),mx2(rip(jj).indstart),'g+');
            plot(tr(rip(jj).indstop),mx2(rip(jj).indstop),'b+');
            
            plot(tr(rip(jj).indallmax),yrx.fz(rip(jj).indallmax)+20,'yO');
            plot(tr(rip(jj).indallmin),yrx.fz(rip(jj).indallmin)+20,'mO');
            
            plot(tr(rip(jj).indabsmax),yrx.fz(rip(jj).indabsmax)+20,'r>');
            plot(tr(rip(jj).indabsmin),yrx.fz(rip(jj).indabsmin)+20,'b<');
            
            
            plot(tr(rip(jj).indstart:rip(jj).indstop),yrx.z(rip(jj).indstart : rip(jj).indstop)+30,'b');
            plot(tr(rip(jj).indstart:rip(jj).indstop),yrx.fz(rip(jj).indstart : rip(jj).indstop)+20,'b');
            
            
            hold off
            disp(' ')
            disp(['final ' num2str([ jj tmp2(jj,:)]) ] )
            disp('-')
            
            pause
            
        end
        
    end
end

% RG - avoids error when no ripples detected
if ~exist('rip','var')
    rip = [];
    allres = [];
    return
end
% /RG



% final bits of analysis---------------------------------------------------
% distribution of values:
allres.lenbin=[0:0.005:0.2];
allres.ampbin=[0:10:300];
allres.freqbin=[120:5:250];

allres.len=histc([rip.length],allres.lenbin);
allres.mlen=mean([rip.length]);
allres.stdlen=std([rip.length]);

allres.freq=histc([rip.freq],allres.freqbin);
allres.mfreq=mean([rip.freq]);
allres.stdfreq=std([rip.freq]);


allres.amp=histc([rip.amplitude],allres.ampbin);
allres.mamp=mean([rip.amplitude]);
allres.stdamp=std([rip.amplitude]);


% ripple props over time (10 sec bins)
%disp('Hist')
allres.tsbin=[tr(1):20:tr(end)];
[allres.c_time,allres.ind_time]=histc([rip.tsstart],allres.tsbin);

for i=1:length(allres.tsbin)
    if allres.c_time(i)>0
        allres.amp_time(i)=mean([rip([allres.ind_time]==i).amplitude]);
        allres.len_time(i)=mean([rip([allres.ind_time]==i).length]);
        allres.freq_time(i)=mean([rip([allres.ind_time]==i).freq]);
    else
        allres.amp_time(i)=0;
        allres.len_time(i)=0;
        allres.freq_time(i)=0;
    end
end
win=5;
[allres.avc] = moving_average(allres.c_time,win,2);
[allres.avamp] = moving_average(allres.amp_time,win,2);
[allres.avlen] = moving_average(allres.len_time,win,2);
[allres.avfreq] = moving_average(allres.freq_time,win,2);


% Make STA of 50 highest amplitude spindles (or all spindles if fewer than
% 50 detected)-RG
allamps=[rip.amplitude];
indallamps=[rip.indabsmax];
[~,ind]=sort(allamps,'descend');

% if length(allamps) >= 50
%     [allres.sta_s,allres.sta_t,allres.sta_E] = sta(cts(indallamps(ind(1:50))) ,cscx.f, cts,'n',[],[],[-1 1],1);
% else
%     [allres.sta_s,allres.sta_t,allres.sta_E] = sta(cts(indallamps),cscx.f, cts,'n',[],[],[-1 1],1);
%     disp('Fewer than 50 spindle events detected!')
% end

% RG
allres.timeranges = opt.ts;

%plot the results

if opt.plt
    fh=figure;
    
    if opt.plt==2, set(fh,'visible','off'); end
    
    set(fh,'name', ' Ripple extraction, summary of results','color','w','PaperPositionMode','auto', 'position',[100 100 600 800])
    subplot(5,2, 1:2)
    plot([0 1],[0 1],'w.');
    set(gca, 'box', 'off', 'xcolor', 'w' , 'ycolor' ,'w','xlim',[0 1])%,...
    %'xtick',[],'ytick',[],'outerposition',proppos);
    xpos=0.01; xpos2=0.5; xpos3=0.8;
    %xpos2=0.5;
    fontnm='FixedWidth';
    text(xpos,1.3, 'Ripples' , 'fontsize',14,'FontName',fontnm)
    text(xpos,0.8, ['Detected: ' num2str(size(rip,2)), ' Ripple events in ' num2str(length(tr)/fsr/60) ' mins' ...
        ' = ' num2str(size(rip,2)/(length(tr)/fsr)) ' Hz ' ], 'fontsize',10,'FontName',fontnm)
    
    text(xpos,0.6, ['mean length: ' num2str(allres.mlen*1000), ' ms' ], 'fontsize',10,'FontName',fontnm)
    text(xpos2,0.6, ['mean amp: ' num2str(allres.mamp), ' uV' ], 'fontsize',10,'FontName',fontnm)
    text(xpos,0.4, ['mean freq: ' num2str(allres.mfreq), ' Hz' ], 'fontsize',10,'FontName',fontnm)
    %text(xpos,0.8, ['      ' fstr.day ' at ' fstr.time ], 'fontsize',14,'FontName',fntstr)
    
    h=subplot(5,2,3:4);
    plot((tr-tr(1))/60,yr,'k');
    
    hold on
    plot((tr([rip.indabsmax])-tr(1))/60,yr([rip.indabsmax]),'*r');
    hold off
    axis tight
    ylabel('uV')
    box off
    set(h,'xticklabel',[],'TickDir','out', 'xcolor','w')
    
    
    
    h=subplot(5,2,5:6);
    plot((allres.tsbin-tr(1))/60,zscore(allres.avc),'k'); hold on
    plot((allres.tsbin-tr(1))/60,zscore(allres.avamp),'b');
    plot((allres.tsbin-tr(1))/60,zscore(allres.avlen),'r');
    plot((allres.tsbin-tr(1))/60,zscore(allres.avfreq),'g');
    hold off
    axis tight
    box off
    legend({'# events', 'Amplitude','Lengths','Frequency'},'FontSize',7,...
        'position',[  0.759  0.761        0.171     0.079])
    xlabel('Time(min)')
    ylabel('z-score')
    set(h,'TickDir','out');
    set(h, 'position',[0.13    0.49  0.775  0.124]);
    
    
    %     subplot(5,2,7)
    %     plot(allres.sta_t,allres.sta_s,'r')
    %     title('Wave triggered average')
    %     axis([-0.05 0.05 -500 500]);
    %     box off
    
    
    subplot(5,2,8)
    line([ allres.mfreq allres.mfreq],[ 0 max(allres.freq)+5],'color',[1 0 0] ,'linewidth',2)
    hold on
    bar(allres.freqbin,allres.freq,'k', 'barwidth',1)
    hold off
    axis([ 120 250 0 max(allres.freq)+5])
    title ('Frequency distribution')
    xlabel('Hz')
    
    
    subplot(5,2,9)
    line([ allres.mamp allres.mamp],[ 0 max(allres.amp)+5],'color',[1 0 0] ,'linewidth',2)
    hold on
    bar(allres.ampbin,allres.amp,'k', 'barwidth',1)
    hold off
    axis([ 0 400 0 max(allres.amp)+5])
    title ('Amplitude distribution')
    xlabel('uV')
    
    
    subplot(5,2,10)
    line([ allres.mlen allres.mlen],[ 0 max(allres.len)+5],'color',[1 0 0] ,'linewidth',2)
    hold on
    bar(allres.lenbin,allres.len,'k','barwidth',1)
    hold off
    axis([ 0 0.2 0 max(allres.len)+5])
    title ('Length distribution')
    
end

end


function [cscx]=rectify_eeg(csc,ts,band,timeranges)

% Modified version of Ulli's function that z-score normalizes using mean
% and std ONLY calculated from the specified time windows

cFs = 1/(ts(2)-ts(1));

if size(csc,1)>1
    csc=csc';
end

cscx.t = ts;
cscx.band = band;
cscx.raw = csc;

% Band-pass filter
cscx.f= eegfilt(csc,cFs,band(1),band(2));

% Rectify
cscx.r= sqrt((cscx.f).^2 );

% Z-score normalize, using only specified time ranges

% Find LFP samples that fall within time windows
valid = rg.helpers.isinrange(ts,timeranges);

sd = std(csc(valid));
mn = mean(csc(valid));

cscx.z = (csc - mn) / sd;

sd = std(cscx.f(valid));
mn = mean(cscx.f(valid));

cscx.fz     = (cscx.f - mn) / sd;
cscx.fzmean = mn;
cscx.fzstd  = sd;

sd = std(cscx.f(valid));
mn = mean(cscx.f(valid));

cscx.rz     = (cscx.r - mn) / sd;
cscx.rzmean = mn;
cscx.rzstd  = sd;

end

function [yR,tR,fsR]=resample_eeg(y,t,fs,fsR)
q = floor(fs/fsR);
yR = resample(y, 1, q);
tR = downsample(t, q);
fsR = fs/q;
end