% MAKE CELL PROFILES
% JULY 23, 2020. 
% Jordan Carpenter

%% add everything to path
warning('off','all') % disable warnings for now..
addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));
addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001"));
addpath(genpath("C:\Users\17145\Documents\github_local\buzcode"));
addpath(genpath("C:\Users\17145\Documents\github_local\mlib6"));
addpath(genpath("C:\Users\17145\Documents\github_local\egoHC"));
%% get spikeAngle for all cells (in degrees)

spikeAngle = cell(1,length(SpikeTrain));
%spikeSpeed = cell(1,length(SpikeTrain));
%spikeAcc = cell(1,length(SpikeTrain));

for sessNum = 1:length(SpikeTrain)
    SHD = cell(1,length(SpikeTrain{1,sessNum}));
    SS = cell(1,length(SpikeTrain{1,sessNum}));
    SA = cell(1,length(SpikeTrain{1,sessNum}));
    for unit = 1:length(SpikeTrain{1,sessNum})
        SHD{1,unit} = getSpikeAngle(hd{1,sessNum}, SpikeTrain{1,sessNum}{1,unit});
        % SS{1,unit} = getSpikeSpeed(speed{1,sessNum}, SpikeTimes{1,sessNum}{1,unit});
        % SA{1,unit} = getSpikeSpeed(accel{1,sessNum}, SpikeTimes{1,sessNum}{1,unit});
    end
    spikeAngle{1,sessNum} = SHD;
    % spikeSpeed{1,sessNum} = SS;
    % spikeAcc{1,sessNum} = SA;
end

% Get home/random well locations for all sessions
[hwLoc, rdLoc] = getWellLoc(labNotes, trialType);






%% Generate cell profile figures

for sessNum = 30%:length(SpikeTrain)
    
    if ~isempty(SpikeTimes{1,sessNum}) % skip empty trials
        
        % Grab some info about current *session*
        sampleRate = mode(diff(pos{1,sessNum}(:,1)));% video-tracking sampling frequency (S)
        date = sessInfo{1,sessNum}.trial_date;
        time = sessInfo{1,sessNum}.trial_time;
        P = pos{1,sessNum};
        
        % Grab window limits for pos tracking
        xMin = str2double(sessInfo{1,sessNum}.window_min_x{1,1})-10;
        xMax = str2double(sessInfo{1,sessNum}.window_max_x{1,1})+10;
        yMin = str2double(sessInfo{1,sessNum}.window_min_y{1,1})-10;
        yMax = str2double(sessInfo{1,sessNum}.window_max_y{1,1})+10;
        
        
        for unit = 1:length(SpikeTrain{1,sessNum})
            
            % Grab some info about current *neuron*
            ST = SpikeTimes{1,sessNum}{1,unit}; 
            UID = UniqueID{1,sessNum}{1,unit}; 
            SPK_A = spikeAngle{1,sessNum}{1,unit};
            % SSPD = spikeSpeed{1,sessNum}{1,unit};
            
            
            %% Calculate some stuff
            [spkPos, spkInd] = data.getSpikePositions(ST,P); % why do i have this?
            map = analyses.map(P, ST, 'smooth', 15, 'binWidth', 4); % FR map
            tc_HD = analyses.turningCurve(spikeAngle{1,sessNum}{1,unit}, P, sampleRate, 'binWidth', 10); %HD tuning curve
            % tc_SPD = analyses.turningCurve(SSPD, P, sampleRate, 'binWidth', 10); %HD tuning curve
            zSpd = [P(:,1), speed{1,sessNum}(:,1)];
            mapSpd = analyses.map(P, zSpd, 'smooth', 15, 'binWidth', 4); % FR map
            [allSpdCnts, hiSpdCnts, loSpdCnts, histEdges] = FRhist(P, ST, speed{1,sessNum});
            
            
            %% Plot everything
            fig = figure('units','normalized','outerposition',[0 0 1 1]); % make fullscreen fig
            set(gcf,'color','w');
            dateStr = char(extractBetween(sessInfo{1,sessNum}.trial_date{1,1},",","20"));
            dateStr = dateStr(find(~isspace(dateStr)));
            timeStr = sessInfo{1,sessNum}.trial_time{1,1}(1:5);
            figTit = strcat('UID:', sprintf('%.f', UID), '\SESS:', sprintf('%.f', sessNum), '\DAT:', dateStr, '\TIM:', timeStr, '\TYPE:', trialType{1,sessNum});
            fig.Name = figTit; % set figure name
            
            % PATHPLOT (STANDARD)
            subplot(5,5,1)
            pathPlot(pos{1,sessNum},SpikeTimes{1,sessNum}{1,unit})
            title("Path Plot")
            xlabel("x")
            ylabel("y")
            xlim([xMin xMax])
            ylim([yMin yMax])
            box off
            
            % PATHPLOT (HD)
            subplot(5,5,2)
            pathPlot_HD(pos{1,sessNum},SpikeTimes{1,sessNum}{1,unit}, hd{1,sessNum})
            colormap(gca,'hsv')
            caxis([0 360])
            colorbar
            xlabel("x")
            ylabel("y")
            xlim([xMin xMax])
            ylim([yMin yMax])
            box off
            
            % FR MAP
            subplot(5,5,3)
            % imagesc(map.z); need to uninvert this if its gonna be used
            plot.colorMap(map.z)
            colorbar
            colormap(gca,'jet')
            % title("Rate Map")
            xlabel("xbins")
            ylabel("ybins")
            box off
            
            % HD PATHPLOT
            subplot(5,5,7)
            scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],hd{1,sessNum}(:,1),'.')
            colormap(gca,'hsv')
            caxis([0 360])
            colorbar
            xlim([xMin xMax])
            ylim([yMin yMax])
            title("HD Path")
            xlabel("x")
            ylabel("y")
            box off         
            
            % ACCEL PATHPLOT
            subplot(5,5,12)
            scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],accel{1,sessNum}(:,1),'.')
            colorbar
            colormap(gca,'copper')
            xlim([xMin xMax])
            ylim([yMin yMax])
            caxis([0 200])
            title("Accel Path")
            box off
            
            % SPEED PATHPLOT
            subplot(5,5,17)
            scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],speed{1,sessNum}(:,1),'.')
            colorbar
            colormap(gca,'copper')
            xlim([xMin xMax])
            ylim([yMin yMax])
            caxis([0 200])
            title("Speed Path")
            box off
            
            
            % HD TUNING CURVE
            subplot(5,5,8)
            nanSum = sum(isnan(hd{1,sessNum}));
            HDnans = sprintf('%.f', nanSum);
            plot(tc_HD(:,1), tc_HD(:,2), 'Color', 'k', 'LineWidth', 1.5)
            t = text(250, .5, strcat('NaNs=', HDnans));
            t.FontSize = 6;
            xlim([0 360])
            title("HD TC")
            xlabel("head angle (deg)")
            box off
            % ylabel("TC value")
            
            % SPEED TUNING CURVE
            subplot(5,5,18)
            [spkSpd] = getSpikeSpeed(P,speed{1,sessNum},ST); % 10 bins
            
            % subplot 8 is dependent on trial type
            if string(trialType{1,sessNum}) == "FM" && isstruct(fmEvents{1,sessNum})
                
                % For FM trials:
                subplot(5,5,9)
                taskPhz = parseTask(fmEvents{1,sessNum}.events, pos{1,sessNum});
                [TC_taskPhz] = taskPhz_tuningCurve(pos{1,sessNum}, SpikeTrain{1,sessNum}{1,unit}, taskPhz);
                for phz=1:length(TC_taskPhz)
                    if isnan(TC_taskPhz(phz))
                        TC_taskPhz(phz) = 0;
                    end
                end
                plot([1 2 3 4 5 6], TC_taskPhz, 'LineWidth', 1.5, 'Color', 'k')
                xticks([1 2 3 4 5 6])
                xticklabels({'SON', 'HW', 'F', 'RW', 'TH', 'SOF'})
                title("Task Phase TC")
                xlabel("task phase")
                ylabel("fr (Hz)")
                box off
                
            elseif string(trialType{1,sessNum}) == "OF"
                
                % For OF trials:
                subplot(5,5,9)
                plot(1:3, 1:3,'LineWidth', 1.5, 'Color', 'k') % plot random thing (for now)
                title("OF SES")
                box off
                
            end
            
            % FR HISTOGRAM- ALL SPEEDS    
            subplot(5,5,6)
            histogram('BinEdges',histEdges,'BinCounts',allSpdCnts, 'FaceColor', 'k')
            title("FR Distrib")
            ylim([0 max(allSpdCnts)])
            xlabel("counts")
            xlabel("fr (Hz)")
            box off
            
            % FR HISTOGRAM- BELOW 5 CM/S
            subplot(5,5,11)
            histogram('BinEdges',histEdges,'BinCounts',loSpdCnts, 'FaceColor', 'k')
            title("FR Distrib: spd<5cm/s")
            ylim([0 max(allSpdCnts)])
            ylabel("counts")
            xlabel("fr (Hz)")
            box off

            
            % FR HISTOGRAM- ABOVE 5 CM/S
            subplot(5,5,16)
            histogram('BinEdges',histEdges,'BinCounts',hiSpdCnts, 'FaceColor', 'k')
            title("FR Distrib: spd>5cm/s")
            ylim([0 max(allSpdCnts)])
            ylabel("counts")
            xlabel("fr (Hz)")
            box off
            
            % PSTH + RASTER (if fmEvents exists)
            
            if string(trialType{1,sessNum}) == "FM" && isstruct(fmEvents{1,sessNum})
                
                % FOR FM TRIALS:
                binsz = 100; % 100 ms
                pre = 5000; % 5000 ms
                post = 5000; % 5000 ms
                [homeEvnts, randEvnts] = trigTimes(fmEvents{1,sessNum}.events);

                
                % FOR HOME DRINKING EVENTS:
                [psth_home trialspx_home] = mpsth(ST, homeEvnts, 'pre', 5000, 'post', 5000, 'binsz', 100);
                [rastmat_home timevec_home] = mraster(trialspx_home,pre,post);

                % PSTH (HOME EVENTS)
                subplot(5,5,5)
                bar(psth_home(:,1)+binsz,psth_home(:,2),'k','BarWidth',1)
                axis([min(psth_home(:,1))-10 max(psth_home(:,1))+10 0 max(psth_home(:,2))+1])
                xlabel('peri-stimulus time'),ylabel(['cntsPer ' num2str(binsz) 'ms bin / fr (Hz)']);
                title("PSTH: Home Events")
                hold on
                xline(0,'--r');
                hold off
                box off
                
                % RASTER (HOME EVENTS)
                subplot(5,5,10)
                for i = 1:numel(trialspx_home)
                    plot(timevec_home,rastmat_home(i,:)*i,'Color','k','Marker','.','MarkerSize',3,'LineStyle','none')
                    hold on
                end
                axis([-pre+10 post+10 0.5 numel(trialspx_home)+0.5])
                xlabel('time (ms)'),ylabel('trials')
                xlabel('peri-stimulus time (ms)')
                title("Raster: Home Events")
                hold on
                xline(0,'--r');
                hold off
                box off
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % FOR RANDOM DRINKING EVENTS:
                [psth_rand trialspx_rand] = mpsth(ST, randEvnts, 'pre', 5000, 'post', 5000, 'binsz', 100);
                [rastmat_rand timevec_rand] = mraster(trialspx_rand,pre,post);
  
                % PSTH (RAND EVENTS)
                subplot(5,5,15)
                bar(psth_rand(:,1)+binsz,psth_rand(:,2),'k','BarWidth',1)
                axis([min(psth_rand(:,1))-10 max(psth_rand(:,1))+10 0 max(psth_rand(:,2))+1])
                xlabel('peri-stimulus time'),ylabel(['cntsPer ' num2str(binsz) 'ms bin / fr (Hz)']);
                title("PSTH: Random Events")
                hold on
                xline(0,'--r');
                hold off
                box off
                
                % RASTER (RAND EVENTS)
                subplot(5,5,20)
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
                
                 % DISTANCE FROM HW TC
                subplot(5,5,14)
                objDist(P, hwLoc{1,sessNum}, ST)
            
                % GOAL DIRECTION TC (HD)- SAREL (??)
                subplot(5,5,19)
                goalDirSar(P, hwLoc{1,sessNum}, hd{1,sessNum}, ST);
                
            elseif string(trialType{1,sessNum}) == "OF" % for open field sessions
                
                for pltNum = [5 10 14 15 19 20]
                    % For OF trials:
                    subplot(5,5,pltNum)
                    plot(1:3, 1:3, 'Color', 'k', 'LineWidth', 1.5) % plot random thing (for now)
                    title("OF SES")
                    box off
                end
                
            else
                for pltNum = [5 10 14 15 19 20]
                    % For TS/other trials:
                    subplot(5,5,pltNum)
                    plot(1:3, 1:3, 'Color', 'k', 'LineWidth', 1.5) % plot random thing (for now)
                    title("OTHER SES")
                    box off
                end
                
            end  
            
            % THETA PHASE TC
            subplot(5,5,4)
            [thetaPhz, spkThetaPhz, tc_TP] = getThetaPhz(rawEEG{1,sessNum}{1,1}, ST);
            plotThetaPhzTC(tc_TP)
            
            % ACCEL TC
            subplot(5,5,13)
            [spkAcc] = getSpikeAccel(P, accel{1,sessNum}, ST);
            
            % click through all figures
            pause
            clf 
        
        end
    end
end

warning('on','all') % re-enable all warning messages for future use