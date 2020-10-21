function [hd_occ, allo_occ, ego_occ, time_occ] = get_binned_occupancy(position, refLoc, whichPlot)
%PLOT_BINNED_OCCUPANCY Summary of this function goes here
%   Detailed explanation goes here


% parse position vector
t = position(:,1);
x1 = position(:,2); y1 = position(:,3);
x2 = position(:,4); y2 = position(:,5);

% compute head direction
head_direction = get_hd(position);

% compute 2D spatial bins
nBins = 10; 
sampleRate = mode(diff(t));

% Compute spatial occupancy and indices for X and Y bins (binX & binY)
[~,xEdges,yEdges,binX,binY] = histcounts2(x1,y1,nBins);

%% calculate egocentric bearing (for all timepoints)
% find midpoint between two LEDs
midX = (x1+x2)/2; midY = (y1+y2)/2;

% break down refLoc
rlX = refLoc(1,1); rlY = refLoc(1,2);

% find allocentric + egocentric 'bearing' at each timepoint (AB/EB)
alloAng = rem(atan2d(rlY-midY, rlX-midX)+180, 360);
egoAng = alloAng - head_direction;
% correct for negative angles (egoAng)
neg_idx = find(egoAng<0);
egoAng(neg_idx) = egoAng(neg_idx)+360;

%% plot

% initalize cell arrays
hd_occ = cell(10,10);
allo_occ = cell(10,10);
ego_occ = cell(10,10);
time_occ = cell(10,10);

hist_bins = 20; % number of bins for histogram

fH = figure('Visible', 'on');
set(gcf,'color','w');

count = 1;
for xx = 1:nBins
    for yy = 1:nBins
        indices = find(xx == binX & yy == binY);
        time_occ{xx,yy} = length(indices)*sampleRate; % occupancy (s)
        
        % calculate values for current 2D spatial bin
        hd_occ{xx,yy} = head_direction(indices); 
        allo_occ{xx,yy} = alloAng(indices); 
        ego_occ{xx,yy} = egoAng(indices);
            
        % choose which measure to plot
        switch whichPlot
            case "hd"
                sDist = hd_occ{xx,yy};
            case "allo"
                sDist = allo_occ{xx,yy};
            case "ego"
                sDist = ego_occ{xx,yy};
        end
        
        % plot
        subplot(10,10,count)
        h{xx,yy} = polarhistogram(sDist,20,'FaceAlpha',.3);
        pax{xx,yy} = gca;
        thetaticks([0]);
        rticks([(nanmax(h{xx,yy}.Values))]);
%         plot_title = strcat('x:', sprintf('%.f', xx), ', y:', sprintf('%.f', yy));
%         title(plot_title);
        
        count = count + 1;
        
    end
end

end

%% SCRATCH CODE
%         subAx{xx,yy} = subplot(10, 10, count, polaraxes);
%         obj{xx,yy} = CircHist(sDist, hist_bins, 'parent', subAx{xx,yy});
%         thetaticks(obj{xx,yy}.polarAxs, 0:90:360);
%         obj{xx,yy}.drawScale; % update scale bar
%         obj{xx,yy}.polarAxs.ThetaZeroLocation = 'right';
%         rl = rlim(obj{xx,yy}.polarAxs);
%         delete(obj{xx,yy}.rH)
%         obj{xx,yy}.drawArrow(obj{xx,yy}.avgAng, obj{xx,yy}.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')
%         plotName = strcat('x=', sprintf('%.f', xx), ',y=', sprintf('%.f', yy));        
%         title(plotName)

% % make the polar plots
% count = 1;
% for xx = 1:nBins
%     for yy = 1:nBins
%     maxVal(count) = nanmax(h{xx,yy}.Values);
%     count = count + 1;
%     end
% end
% 
% % find maximum value
% maxmax = nanmax(maxVal);
% 
% % change max value for all plots (normalize)
% count = 1;
% for xx = 1:nBins
%     for yy = 1:nBins
%     pax{xx,yy}.RLim = [1, maxmax];
%     count = count + 1;
%     end
% 


% %% compute spatial bin centers 
% for i = 1:length(xEdges)
%     if i+1 <= length(xEdges)
%         xCenter(i) = ((xEdges(i+1)-xEdges(i))/2)+xEdges(i);
%     end
% end
% 
% for i = 1:length(yEdges)
%     if i+1 <= length(yEdges)
%         yCenter(i) = ((yEdges(i+1)-yEdges(i))/2)+yEdges(i);
%     end
% end
% 
% % make a vector of the bin centers
% count = 1;
% for xx = 1:length(xCenter)
%     for yy = 1:length(yCenter)
%         ctrLocs(count,1:2) = [xCenter(xx), yCenter(yy)];
%         count = count+1;
%     end
% end
% 
% % reshape the bin centers
% binCtrs=[];
% number = 1;
% for ii = 1:10
%     binCtrs= [binCtrs; binCenters(ii:10:end,:)];
%     number = number + 1;
% end



