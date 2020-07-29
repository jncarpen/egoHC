% MAKE CELL PROFILES
% JULY 23, 2020.

%% get spikeAngle for all cells (in degrees)

spikeAngle = cell(1,length(SpikeTrain));
spikeSpeed = cell(1,length(SpikeTrain));
spikeAcc = cell(1,length(SpikeTrain));

for sessNum = 1:length(SpikeTrain)
    SHD = cell(1,length(SpikeTrain{1,sessNum}));
    SS = cell(1,length(SpikeTrain{1,sessNum}));
    SA = cell(1,length(SpikeTrain{1,sessNum}));
    for unit = 1:length(SpikeTrain{1,sessNum})
        SHD{1,unit} = getSpikeAngle(hd{1,sessNum}, SpikeTrain{1,sessNum}{1,unit});
        SS{1,unit} = getSpikeSpeed(speed{1,sessNum}, SpikeTrain{1,sessNum}{1,unit});
        SA{1,unit} = getSpikeSpeed(accel{1,sessNum}, SpikeTrain{1,sessNum}{1,unit});
    end
    spikeAngle{1,sessNum} = SHD;
    spikeSpeed{1,sessNum} = SS;
    spikeAcc{1,sessNum} = SA;
end



%% trialType (still empty stuff)

% trialType = cell(1,length(sessType));
% for sessNum = 1:length(sessType)
%     if ~isempty(sessType{1,sessNum})
%         binaryArr = cellfun(@isempty,test); % empty cells represented as 0.
%         idx = find(binaryArr == 0, 1); % first instance of a non-empty array
%         trialType{1,sessNum} = sessType{1,sessNum}{1,idx};
%     end
% end



%% Generate cell profile figures

for sessNum = 28%:length(SpikeTrain)
    
    if ~isempty(SpikeTrain{1,sessNum}) % skip empty trials
        
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
        
        % Bin SPEED
%         customCM = zeros(length(speed{1,sessNum}),3);
%         speedBins = 20;
%         speedEdges = linspace(0,25,speedBins+1);
%         [speedOccupancy,edgesHist,speedInds] = histcounts(speed{1,sessNum}(:,1), speedEdges);
%         cmap = jet(speedBins+1);
%         for bin = 0:speedBins
%             SpdIdx = find(speedInds == bin);
%             customCM(SpdIdx,1:3) = cmap(bin+1);
%         end
        
        
        
        
        for unit = 1:length(SpikeTrain{1,sessNum})
            
            % Grab some info about current *neuron*
            ST = SpikeTimes{1,sessNum}{1,unit}; 
            UID = UniqueID{1,sessNum}{1,unit}; 
            SPK_A = spikeAngle{1,sessNum}{1,unit};
            SSPD = spikeSpeed{1,sessNum}{1,unit};
            
            
            % Calculate some stuff
            % [spkPos, spkInd] = data.getSpikePositions(ST,P); % why do i have this?
            map = analyses.map(P, ST, 'smooth', 15, 'binWidth', 4); % FR map
            tc_HD = analyses.turningCurve(SPK_A, P, sampleRate, 'binWidth', 10); %HD tuning curve
            tc_SPD = analyses.turningCurve(SSPD, P, sampleRate, 'binWidth', 10); %HD tuning curve
            
            zSpd = [P(:,1), speed{1,sessNum}(:,1)];
            mapSpd = analyses.map(P, zSpd, 'smooth', 15, 'binWidth', 4); % FR map

            
            % Plot everything
            figure
            
            
            % PATHPLOT (STANDARD)
            subplot(3,3,1)
            pathPlot(pos{1,sessNum},SpikeTimes{1,sessNum}{1,unit})
            title("Path Plot")
            xlabel("X")
            ylabel("Y")
            xlim([xMin xMax])
            ylim([yMin yMax])
            box off
            
            % PATHPLOT (HD)
            subplot(3,3,2)
            pathPlot_HD(pos{1,sessNum},SpikeTimes{1,sessNum}{1,unit}, hd{1,sessNum})
            xlabel("X")
            ylabel("Y")
            xlim([xMin xMax])
            ylim([yMin yMax])
            box off
            
            
            % FR MAP
            subplot(3,3,3)
            % imagesc(map.z); need to uninvert this if its gonna be used
            plot.colorMap(map.z)
            colorbar
            colormap(gca,'jet')
            % title("Rate Map")
            xlabel("xBins")
            ylabel("yBins")
            box off
            
            
            % ACCEL PATHPLOT
            subplot(3,3,4)
            scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],accel{1,sessNum}(:,1),'.')
            colorbar
            colormap(gca,'summer')
            xlim([xMin xMax])
            ylim([yMin yMax])
            caxis([0 200])
            title("Accel Path")
            box off

            
            % SPEED PATHPLOT
            subplot(3,3,5)
            scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],speed{1,sessNum}(:,1),'.')
            colorbar
            colormap(gca,'summer')
            xlim([xMin xMax])
            ylim([yMin yMax])
            caxis([0 200])
            title("Speed Path")
            box off
            
            % HD PATHPLOT
            subplot(3,3,6)
            scatter(pos{1,sessNum}(:,2), pos{1,sessNum}(:,3),[10],hd{1,sessNum}(:,1),'.')
            colorbar
            colormap(gca,'hsv')
            xlim([xMin xMax])
            ylim([yMin yMax])
            caxis([0 360])
            title("HD Path")
            box off
            
            
            % HD TUNING CURVE
            subplot(3,3,7)
            plot(tc_HD(:,1), tc_HD(:,2), 'Color', 'k', 'LineWidth', 1.5)
            xlim([0 360])
            title("HD Tuning Curve")
            xlabel("Head Angle (deg)")
            box off
            % ylabel("TC value")
            
            % SPEED TUNING CURVE
            subplot(3,3,8)
            plot(tc_SPD(:,1), tc_SPD(:,2), 'Color', 'k', 'LineWidth', 1.5)
            title("Speed Tuning Curve")
            xlabel("Speed")
            xlim([0 200])
            box off
            % ylabel("TC value")
            % subplot 8 is dependent on trial type
            if string(trialType{1,sessNum}) == "FM" && isstruct(fmEvents{1,sessNum})
                
                % For FM trials:
                subplot(3,3,9)
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
                title("Task Phase")
                box off
                
            elseif string(trialType{1,sessNum}) == "OF"
                
                % For OF trials:
                subplot(3,3,9)
                plot(1:3, 1:3)
                title("OF SES")
                box off
                
            end
            
        end
    end
end