
% From: Oyvind & Sebastian
% I will modify this and try to use it on my data
% JC

function [Vectormap,r_bin,n,orientationCurve,polarbins,circlebins] = makeVecMaps(pos_,SpikeTimes_,objectPos_)

%   INPUTS
%   'pos':                  [t x1 y1 x2 y2]
%   'spike2posInd':         positions of animal when cell spikes
%   'objectPos':            [x y]
%
%   OUTPUTS
%   'Vectormap':
%   'r_bin':
%   'n':
%   'orientationCurve':
%   'polarbins':
%   'circlebins':

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define parameters
n = 36; % number of angular bins
r_bin = 3; %5 from beginning
smoothF = 1.5; %1.5;
maxDistMap = 12; % max distance threshold?
%startwidth = 1;

% get pos indices for spike times
spike2posInd = knnsearch(pos_(:,1),SpikeTimes_);

 
%% find distance and direction of each pos-sample relative to each sample of moving object position.
% Column 1 in polarpos will be orientation of mouse relative to
% object. Column 2 will be distance between object and mouse

for i = 1:size(pos_,1)
  [polarpos(i,1),polarpos(i,2)] = cart2pol(pos_(i,2)-objectPos_(1,1), pos_(i,3)-objectPos_(1,2));                 %[polarpos(i,1),polarpos(i,2)] = cart2pol(pos(i,2)-objectpos(i,2), pos(i,3)-objectpos(i,3));
end

polarpos(:,1) = polarpos(:,1)+pi;

%%
polarbins = linspace(0,2*pi,n+1)';
radius = max(polarpos(:,2));
circlebins = 1:r_bin:radius+r_bin;

for i = 1:length(polarbins)-1
    sectorbins{i} = find(polarpos(:,1)>polarbins(i) & polarpos(:,1)<polarbins(i+1));
end

%for i = 1:length(spike2posInd)
    for j = 1:length(sectorbins)
    oriSpikes{1,j} = intersect(sectorbins{j},spike2posInd);     %changed from curly to normal brackets
    end
%end

%% make tuningcurve for orientation alone

%for i = 1:size(oriSpikes,1)
    for j = 1:size(oriSpikes,2)
        orientationCurve(1,j) = numel(oriSpikes{1,j})./(numel(sectorbins{j})/120);
    end
    orientationCurve(1,:) = general.smooth(orientationCurve(1,:),smoothF);
%end
%%
for i = 1:length(circlebins)-1
    distancebins{i} = find(polarpos(:,2)>circlebins(i) & polarpos(:,2)<circlebins(i+1));
end

%for i = 1:length(spike2posInd)
    for j = 1:length(distancebins)
    distanceSpikes{1,j} = intersect(distancebins{j},spike2posInd);   %changed from curly to normal brackets
    end
%end


%note: some linear indices in spike2posInd may be repeated; intersect will
%remove repeating elements; so it may look like we have "lost" spikes


%%

for i = 1:length(sectorbins)
    for j = 1:length(distancebins)
        posmatrix(i,j) = numel(intersect(cell2mat(sectorbins(i)),cell2mat(distancebins(j))))/120;
    end
end
posmatrix(posmatrix==0) = NaN;

for k = 1:size(oriSpikes,1)
    for i = 1:size(oriSpikes,2)
        for j = 1:size(distanceSpikes,2)
            spikematrix(i,j) = numel(intersect(cell2mat(oriSpikes(k,i)),cell2mat(distanceSpikes(k,j))));
        end
    end
   % Smooth vectormaps with Gaussian kernel;
   C = spikematrix./posmatrix;
   nancurve = find(isnan(C(:)));
   C(isnan(C)) = 0;
   C = general.smooth(C,smoothF);
   sizeC = size(C);
   C = C(:);
   C(nancurve)= NaN;
   C = reshape(C,sizeC);
   measuredDist = size(C,2);
   if measuredDist<maxDistMap
       C(:,measuredDist+1:maxDistMap) = NaN;
   end
   C = C(:,1:maxDistMap);

   spikeMatrix{k} = spikematrix; 
   Vectormap = C;                             %Vectormap{k} = C;
end



% %% Plot
% 
% 
% 
figure();
h = pcolor(C);
%remove black edges from figure (Matlab standard)
set(h, 'EdgeColor', 'none');

axis square
axis([1 size(C,2) 1 size(C,1)])

%format and label x- and y axes

%distance goes from 0 to 60 cm (x-axis)
xticks = [1 size(C,2)/2 size(C,2)];
xticklabels = [0 round(r_bin*xticks(2)) round(r_bin*xticks(3))];
set(gca,'XTick',xticks,'XTickLabel', xticklabels) 
xlabel('Distance (cm)');

%orientation goes from 0 to 360 deg (y-axis)
yticks = [1 size(C,1)/2 size(C,1)];
yticklabels = [0 180 360];
set(gca,'YTick',yticks,'YTickLabel', yticklabels);
ylabel('Orientation (deg)'); 


%check with Øyvind what this one was 
%plot([0 size(C,2)],[mloc mloc],'w--','linewidth',1)   




% k = 1;
% while k <= length(spikeMatrix)
% figure
% C = Vectormap{1,k};
% h = pcolor(C);
% set(h, 'EdgeColor', 'none');
% %title([dateSession, ' Cell ',num2str(k)])
% axis square
% axis([1 size(C,2) 1 size(C,1)])
% plot([0 size(C,2)],[mloc mloc],'w--','linewidth',1)
% yticks = [1 size(C,1)/2 size(C,1)];
% yticklabels = [0 180 360];
% xticks = [1 size(C,2)/2 size(C,2)];
% xticklabels = [0 round(r_bin*xticks(2)) round(r_bin*xticks(3))];
% set(gca,'XTick',xticks,'XTickLabel', xticklabels) 
% set(gca,'YTick',yticks,'YTickLabel', yticklabels)
% ylabel('Orientation (deg)')
% xlabel('Distance (cm)')
% 
% figure
% n = pcolor(objectmaps{1, k}.z);
% set(n, 'EdgeColor','None')
% axis square
% 
% k = k + 1;
% end
% 
% figure
% imagesc(ratemaps{1,5}.z)
% axis square









%%

