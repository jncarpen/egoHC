function [hd_occ, allo_occ, ego_occ, time_occ, CVM_Dist] = get_binned_occupancy(position, refLoc, whichPlot)
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

fH = figure('Visible', 'on');
set(gcf,'color','w');

% number of bins for histogram
num_bins = 20;
        
% generate a uniform distribution between 0 and 360
% figure
% a = 0; b = 360;
% uniDist = (b-a).*rand(100000,1) + a;
% uniHist = polarhistogram(deg2rad(uniDist),num_bins,'FaceAlpha',.3);  % output in radians
% data = uniHist.Data;

% % find the average value in each bin
% edges = uniHist.BinEdges'; data = uniHist.Data;
% for bb = 1:length(edges)-1
%     for dd = 1:length(data)
%         data_binned = data(data>edges(bb)&data<=edges(bb+1));
%         uniMean(bb) = nanmean(data_binned);
%     end
% end


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
        
        % find PDF of sDist
        P1 = [];
        sDist_rnd = round(sDist);
        for val = 1:360
            if or(val == 0, val == 360)
                sum_wrap = sum(sDist_rnd == 0) + sum(sDist_rnd == 360);
                P1(360) = sum_wrap/length(sDist_rnd);
            else
                P1(val) = sum(sDist_rnd == val)/length(sDist_rnd);
            end
        end
            
        % define pdf vectors
        X = linspace(1, 360, 360)';
        P1 = P1' + eps;
        P2 = ones(360,1)/360; % PDF of uniform distribution
        
%         dist(xx,yy) = abs(sum(P1-P2))*100000000000000;
        
        % find Kullback-Leibler divergence
%         KL(xx,yy) = kldiv(X,P1,P2);
        CVM_Dist(xx,yy) = Cramer_Von_Mises(P1,P2);
%         WS_Dist(xx,yy) = Wasserstein_Dist(P1,P2);
        
        
        % find distance from uniform distribution
%         dist_from_uni(xx,yy) = abs(deg2rad(mean(sDist))-mean(data));
        

        % plot
        subplot(10,10,count)
        h{xx,yy} = polarhistogram(sDist,num_bins,'FaceAlpha',.3);
        pax{xx,yy} = gca;
        thetaticks([0]);
        rticks([(nanmax(h{xx,yy}.Values))]);
%         uniform_count = repmat(round(sum(h{xx,yy}.Values)/num_bins), 1, num_bins);
%         dist_from_uni(xx,yy) = sum(h{1,1}.Values-uniform_count);
        plot_title = strcat('CVM', sprintf('%.2f', CVM_Dist(xx,yy)));
        title(plot_title);
        
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



