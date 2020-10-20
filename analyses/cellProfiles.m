%% MAKE CELL PROFILES
% DESCRIPTION:
% Create a bunch of plots for each cell!

% Jordan Carpenter (jordannc@stud.ntnu.no)
% Created on: July 23, 2020. 
% Last modified: October 9, 2020.

%% add everything to path
warning('off','all') % disable warnings for now..
addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));
addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));
addpath(genpath("C:\Users\17145\Documents\github_local\buzcode"));
addpath(genpath("C:\Users\17145\Documents\github_local\mlib6"));
addpath(genpath("C:\Users\17145\Documents\github_local\egoHC"));


%% get spikeAngle for all cells (in degrees)

% generate SpikeTrain_thresh
[SpikeTrain_thresh] = spkTrain_thresh(pos_cm, SpikeTimes_thresh, SpikeTrain);

spikeAngle = cell(1,length(SpikeTrain));

for sessNum = 1:length(SpikeTrain) % 7 is when training stops
    SHD = cell(1,length(SpikeTrain{1,sessNum}));
    SS = cell(1,length(SpikeTrain{1,sessNum}));
    SA = cell(1,length(SpikeTrain{1,sessNum}));
    for unit = 1:length(SpikeTrain{1,sessNum})
        SHD{1,unit} = getSpikeAngle(hd{1,sessNum}, SpikeTrain_thresh{1,sessNum}{1,unit});

    end
    spikeAngle{1,sessNum} = SHD;
end

% Get home/random well locations for all sessions
[hwLoc, rdLoc] = getWellLoc(labNotes, trialType);


%% Generate necessary values
% this uses SpikeTimes_thresh

for sessNum = 1:length(pos_cm)
    if ~isempty(pos{1,sessNum}) && ~isempty(SpikeTimes_thresh{1,sessNum})
        sessionArray = cell(1,length(SpikeTimes{1,sessNum}));
        for unitNum = 1:length(SpikeTrain{1,sessNum})
            ST_now = SpikeTimes_thresh{1,sessNum}{1,unitNum}; % spifigureke times now
            SpikeTrain_thresh{1,sessNum} = binSpikes(pos_cm{1,sessNum}(:,1), SpikeTimes_thresh{1,sessNum});
        end
    end
end

% get accel_cm
[~, accel_cm] = fix_speed_cm(pos_cm);


%% Generate cell profile figures

for sessNum = 27%1:length(SpikeTrain)
    
    if ~isempty(SpikeTimes{1,sessNum}) % skip empty trials
        disp(sessNum)
        
        % Grab some info about current *session*
        sampleRate = mode(diff(pos{1,sessNum}(:,1)));% video-tracking sampling frequency (S)
        date = sessInfo{1,sessNum}.trial_date;
        time = sessInfo{1,sessNum}.trial_time;
        P = pos_cm{1,sessNum};
        
        % Grab window limits for pos tracking (for pos_cm)
        expansionFactor = 5;
        xMin = nanmin(P(:,2))-expansionFactor;
        xMax = nanmax(P(:,2))+expansionFactor;
        yMin = nanmin(P(:,3))-expansionFactor;
        yMax = nanmax(P(:,2))+expansionFactor;
        
        % get box center coordinates
        ctrCoord = boxCtr{1,sessNum};
        
        % get ref coordinates (home well for FM or center for other)
        refCoord = hwCoord{1,sessNum};
        
        for unit = 4%1:length(SpikeTrain{1,sessNum})
            if length(SpikeTimes_thresh{1,sessNum}{1,unit}) > 80
                disp(unit)
                % Grab some info about current *neuron*
                UID = UniqueID{1,sessNum}{1,unit}; 
                SPK_A = spikeAngle{1,sessNum}{1,unit};
                ST = SpikeTimes_thresh{1,sessNum}{unit};

                %% Calculate tuning curves
                [spkPos, spkInd] = data.getSpikePositions(ST,P); % why do i have this?


                %% Plot everything
                fig = figure('units','normalized','outerposition',[0 0 1 1]); % make fullscreen fig
                set(gcf,'color','w');
                dateStr = char(extractBetween(sessInfo{1,sessNum}.trial_date{1,1},",","20"));
                dateStr = dateStr(find(~isspace(dateStr)));
                timeStr = sessInfo{1,sessNum}.trial_time{1,1}(1:5);
                figTit = strcat('UID:', sprintf('%.f', UID), '\SESS:', sprintf('%.f', sessNum), '\DAT:', dateStr, '\TIME:', timeStr, '\TYPE:', trialType{1,sessNum});
                fileBody = strcat('UID', sprintf('%.f', UID), '_SESS', sprintf('%.f', sessNum), '_TYPE_', trialType{1,sessNum});
                fig.Name = figTit; % set figure name

                % PATHPLOT (STANDARD)
                subplot(4,5,1)
                pathPlot(P, ST)
                pbaspect([1 1 1])
                title("Path Plot")
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                xlim([xMin, xMax])
                ylim([yMin, yMax])
                axis off
                % box off

                % PATHPLOT (HD)
                subplot(4,5,2)
                pathPlot_hd(P, ST, hd{1,sessNum})
                pbaspect([1 1 1])
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                xlim([xMin, xMax])
                ylim([yMin, yMax])
                colormap(gca,'hsv')
                caxis([0 360])
                colorbar
                axis off

                % FIRING RATE MAP
                subplot(4,5,3)
                map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 4); % calculate tuning curve
                peakRate = nanmax(nanmax(map.z));
                rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
                % imagesc(map.z); need to uninvert this if its gonna be used
                plot.colorMap(map.z)
                pbaspect([1 1 1])
                colorbar
                colormap(gca,'jet')
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                title(rate_map_title)
                box off      

                % EGO DISTANCE
                subplot(4,5,7)
                egoDistance(pos_cm, ST, ctrCoord, refCoord, hd, sessNum)

                % SPATIAL OCCUPANCY
                subplot(4,5,11)
                plot.colorMap(map.time, 'ydir', 'normal')
                colormap(gca,'jet')
                pbaspect([1 1 1])
                colorbar
                colormap(gca,'jet')
                title("Occupancy")
                axis off       

                % ALLOCENTRIC PLOT
                subplot(4,5,4)
                [Vectormap,r_bin,n,orientationCurve,polarbins,circlebins] = makeVecMaps(P,ST,refCoord);
                colormap(gca,'jet')
                
                % EGO BEARING
                % LEGEND: blackLine = ctrCoord; redLine = refCoord (hwLoc)
                subplot(4,5,17)
                [tcVals_egoAng] = egoBearing(pos_cm, ST, ctrCoord, refCoord, hd, sessNum, "True", "deg");
               

                % HD TUNING CURVE
                % generate tuning curve values
                SpkAngNow = spikeAngle{1,sessNum}{1,unit}; % get head angle each time neuron spikes
                tc_HD = analyses.turningCurve(SpkAngNow, P, sampleRate, 'smooth', 5, 'binWidth', 10); %HD tuning curve
                tc_HD_even = analyses.turningCurve(SpkAngNow(2:2:end), P(2:2:end), sampleRate, 'smooth', 5, 'binWidth', 10); % even tuning curve
                tc_HD_odd = analyses.turningCurve(SpkAngNow(1:2:end), P(1:2:end), sampleRate, 'smooth', 5, 'binWidth', 10); % even tuning curve
                % plot tuning curves
                subplot(4,5,8)
                plot(tc_HD(:,1), tc_HD(:,2), 'Color', 'r', 'LineWidth', 1.10)
                hold on
                plot(tc_HD_even(:,1), tc_HD_even(:,2), 'Color', 'k', 'LineStyle',':', 'LineWidth', 1.10)
                plot(tc_HD_odd(:,1), tc_HD_odd(:,2), 'Color', 'k', 'LineWidth', 1.10)
                xlim([0 360])
                title("HD TC")
                xlabel("head angle (deg)")
                box off


                % SPEED TUNING CURVE
                subplot(4,5,18)
                [~] = getSpikeSpeed(P, ST);
                
                % subplot 8 is dependent on trial type
                if string(trialType{1,sessNum}) == "FM" && isstruct(fmEvents{1,sessNum})

                    % For FM trials:
                    subplot(4,5,19)
                    taskPhz = parseTask(fmEvents{1,sessNum}.events, P);
                    [TC_taskPhz] = taskPhz_tuningCurve(P, ST, taskPhz);
                    for phz=1:length(TC_taskPhz)
                        if isnan(TC_taskPhz(phz))
                            TC_taskPhz(phz) = 0;
                        end
                    end
                    plot([1 2 3 4 5 6], TC_taskPhz, 'LineWidth', 1.10, 'Color', 'k')
                    xticks([1 2 3 4 5 6])
                    xticklabels({'SON', 'HW', 'F', 'RW', 'TH', 'SOF'})
                    title("Task Phase TC")
                    xlabel("task phase")
                    ylabel("fr (Hz)")
                    box off

                else

                    % For OF trials:
                    subplot(4,5,9)
                    plot(1:3, 1:3,'LineWidth', 1.10, 'Color', 'k') % plot random thing (for now)
                    title("OF SES")
                    box off

                end

                % FR HISTOGRAM- ALL SPEEDS    
                subplot(4,5,16)
                [allSpdCnts, hiSpdCnts, loSpdCnts, histEdges] = FRhist(P, ST, speed_cm{1,sessNum});
                histogram('BinEdges',histEdges,'BinCounts',allSpdCnts, 'FaceColor', 'k')
                title("FR Distrib")
                ylim([0 max(allSpdCnts)])
                xlabel("counts")
                xlabel("fr (Hz)")
                box off
                
                % Spikeplot (heading direction); method 1
                subplot(4,5,6)
                pathPlot_heading(P, ST);
                
                
                % PSTH + RASTER
                if string(trialType{1,sessNum}) == "FM" && isstruct(fmEvents{1,sessNum}) % if the trial is FM

                    % FOR FM TRIALS:
                    binsz = 200; % 100 ms
                    pre = 10000; % 5000 ms
                    post = 10000; % 5000 ms
                    [homeEvnts, randEvnts] = trigTimes(fmEvents{1,sessNum}.events);


                    % FOR HOME DRINKING EVENTS:
                    [psth_home, trialspx_home] = mpsth(ST, homeEvnts, 'pre', pre, 'post', post, 'binsz', 200);
                    [rastmat_home, timevec_home] = mraster(trialspx_home,pre,post);

                    % PSTH (HOME EVENTS)
                    subplot(4,5,5)
                    bar(psth_home(:,1)+binsz,psth_home(:,2),'k','BarWidth',1)
                    axis([min(psth_home(:,1))-10 max(psth_home(:,1))+10 0 max(psth_home(:,2))+1])
                    xlabel('peri-stimulus time'),ylabel(['cntsPer ' num2str(binsz) 'ms bin / fr (Hz)']);
                    title("PSTH: Home Events")
                    hold on
                    xline(0,'--r');
                    hold off
                    box off

                    % RASTER (HOME EVENTS)
                    subplot(4,5,10)
                    for i = 1:numel(trialspx_home)
                        plot(timevec_home,rastmat_home(i,:)*i,'Color','k','Marker','.','MarkerSize',3,'LineStyle','none')
                        hold on
                    end
                    axis([-pre+10 post+10 0.5 numel(trialspx_home)+0.5])
                    xlabel('time (ms)'),ylabel('trials')
                    xlabel('peri-stim time (ms)')
                    title("Raster: Home Events")
                    hold on
                    xline(0,'--r');
                    hold off
                    box off

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % FOR RANDOM DRINKING EVENTS:
                    [psth_rand, trialspx_rand] = mpsth(ST, randEvnts, 'pre', pre, 'post', post, 'binsz', 200);
                    [rastmat_rand, timevec_rand] = mraster(trialspx_rand,pre,post);
                    % PSTH (RAND EVENTS)
                    subplot(4,5,15)
                    bar(psth_rand(:,1)+binsz,psth_rand(:,2),'k','BarWidth',1)
                    axis([min(psth_rand(:,1))-10 max(psth_rand(:,1))+10 0 max(psth_rand(:,2))+1])
                    xlabel('peri-stim time (ms)'),ylabel(['cntsPer ' num2str(binsz) 'ms bin / fr (Hz)']);
                    title("PSTH: Random Events")
                    hold on
                    xline(0,'--r');
                    hold off
                    box off

                    % RASTER (RAND EVENTS)
                    subplot(4,5,20)
                    for i = 1:numel(trialspx_rand)
                        plot(timevec_rand,rastmat_rand(i,:)*i,'Color','k','Marker','.','MarkerSize',3,'LineStyle','none')
                        hold on
                    end
                    axis([-pre+10 post+10 0.5 numel(trialspx_rand)+0.5])
                    xlabel('time (ms)'),ylabel('trials')
                    xlabel('peri-stimulus time (ms)')
                    title("Raster: Random Events")
                    hold on
                    xline(0,'--r');
                    hold off
                    box off

                    % Egocentric bearing TC (wrt heading direction)
                    subplot(4,5,14)
                    % moving_direction_tc(P, ST)
                    [~] = egoBearing_heading(P, ST, ctrCoord, refCoord, "True", "deg", "1");

                    % GOAL DIRECTION TC (HD)- SAREL (??)
                    subplot(4,5,9)
                    goalDirSar_bins = 20;   
                    goalDirSar(P, refCoord, hd{1,sessNum}, ST, goalDirSar_bins); % ref coord is home location for FM and center of box for OF

                elseif string(trialType{1,sessNum}) == "OF" % for open field sessions
                    
                    for pltNum = [5 10 15 19 20]
                        % For OF trials:
                        subplot(4,5,pltNum)
                        plot(1:3, 1:3, 'Color', 'k', 'LineWidth', 1.10) % plot random thing (for now)
                        title("OF SES")
                        box off
                    end

                else % For TS/other trials:
                    for pltNum = [5 10 15 19 20]
                        subplot(4,5,pltNum)
                        plot(1:3, 1:3, 'Color', 'k', 'LineWidth', 1.10) % plot random thing (for now)
                        title("OTHER SES")
                        box off
                    end

                end  

                % THETA PHASE TC
                subplot(4,5,12)
                [thetaPhz, spkThetaPhz, tc_TP] = getThetaPhz(rawEEG{1,sessNum}{1,1}, ST);
                plotThetaPhzTC(tc_TP)

                % ACCEL TC
                subplot(4,5,13)
                [spkAcc] = getSpikeAccel(P, accel_cm{1,sessNum}, ST);

                % save figures
                filename = strcat('D:\egoAnalysis\cell_profiles_oct8\', fileBody, '.png');
                saveas(fig, filename);

                % click through all figures
                % pause
                % close all 
                
            else
                
              % display a message if the neuron does not have enough spikes
              % to pass the criteria
              message = strcat("unit ", sprintf('%.0f', unit), " has too few spikes to include");
              disp(message)
              
            end
        end
    end
end
warning('on','all') % re-enable all warning messages for future use




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCRATCH CODE

%             
%             % FR HISTOGRAM- BELOW 5 CM/S
%             subplot(4,5,11)
%             histogram('BinEdges',histEdges,'BinCounts',loSpdCnts, 'FaceColor', 'k')
%             title("FR Distrib: spd<5cm/s")
%             ylim([0 max(allSpdCnts)])
%             ylabel("counts")
%             xlabel("fr (Hz)")
%             box off


%             % FR HISTOGRAM- ABOVE 5 CM/S
%             subplot(4,5,16)
%             histogram('BinEdges',histEdges,'BinCounts',hiSpdCnts, 'FaceColor', 'k')
%             title("FR Distrib: spd>5cm/s")
%             ylim([0 max(allSpdCnts)])
%             ylabel("counts")
%             xlabel("fr (Hz)")
%             box off


%             % SPEED PATHPLOT
%             subplot(4,5,17)
%             scatter(P(:,2), P(:,3),[10],speed{1,sessNum}(:,1)/2.6,'.')
%             pbaspect([1 1 1])
%             colorbar
%             colormap(gca,'copper')
%             xlim([xMin xMax])
%             ylim([yMin yMax])
%             caxis([0 200])
%             title("Speed Path")
%             box off


%             % HD PATHPLOT
%             subplot(4,5,7)
%             scatter(P(:,2), P(:,3),[10],hd{1,sessNum}(:,1),'.')
%             pbaspect([1 1 1])
%             colormap(gca,'hsv')
%             caxis([0 360])
%             colorbar
%             xlim([xMin xMax])
%             ylim([yMin yMax])
%             title("HD Path")
%             xlabel("x")
%             ylabel("y")
%             box off   


%             % ACCEL PATHPLOT
%             subplot(4,5,12)
%             scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],accel{1,sessNum}(:,1),'.')
%             colorbar
%             colormap(gca,'copper')
%             xlim([xMin xMax])
%             ylim([yMin yMax])
%             caxis([0 200])
%             title("Accel Path")
%             box off
    
% %         get objectPos (depending on session type)
%         if string(trialType{1,sessNum}) == "FM" && hwLoc{1,sessNum} == "36"
%         elseif string(trialType{1,sessNum}) == "FM" && hwLoc{1,sessNum} == "37"
%         else
%         end
        
        % Grab window limits for pos tracking
%         xMin = str2double(sessInfo{1,sessNum}.window_min_x{1,1})-10;
%         xMax = str2double(sessInfo{1,sessNum}.window_max_x{1,1})+10;
%         yMin = str2double(sessInfo{1,sessNum}.window_min_y{1,1})-10;
%         yMax = str2double(sessInfo{1,sessNum}.window_max_y{1,1})+10;


% %  add text to the figure about the number of NaNs
% t = text(250, .5, strcat('NaNs=', HDnans));
% t.FontSize = 6;
% nanSum = sum(isnan(hd{1,sessNum}));
% HDnans = sprintf('%.f', nanSum);

%     put the legend outside of the main figure later
%     legend('all', 'even', 'odd', 'Location','southwest', 'FontSize',5);
%     legend('boxoff')

% i think this is for the speed path plot
% zSpd = [P(:,1), speed_cm{1,sessNum}(:,1)];
% mapSpd = analyses.map(P, zSpd, 'smooth', 15, 'binWidth', 4); % FR map

% % For object distance tuning curve (retired)
% objDist_bins = 20;
% objDist(P, refCoord, ST, objDist_bins)
