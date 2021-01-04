function plot_modelDynamics(P, ST, model, ref_point)

%% speed threshold 
% parse position data
t = P(:,1); % time (seconds)
x = P(:,2); 
y = P(:,3);
x2 = P(:,4);
y2 = P(:,5);

% tSpk = ST;
% position = P;

% speed
[s, ~] = get_speed(P);
s = s(:,1); % grab first column

% speed threshold before we make the spiketrain
% Get speed at time of spike and put into vector SpikeSpeed
SpikeSpeed = interp1 (t, s, ST); %in cm/s

% Set threshold
thr_d= 4; % this is the threshold set in jercog et al. (diff for dMan)       
thr_u= 100;

% Apply threshold 
a=find(SpikeSpeed>thr_d);
b=find(SpikeSpeed<thr_u);

% make position samples NaN
x(find(s<thr_d))=NaN; x(find(s>thr_u))=NaN;
y(find(s<thr_d))=NaN; y(find(s>thr_u))=NaN;
x2(find(s<thr_d))=NaN; x2(find(s>thr_u))=NaN;
y2(find(s<thr_d))=NaN; y2(find(s>thr_u))=NaN;
position = [t, x, y, x2, y2];

% Combined threshold 
c=intersect(a,b);

% Vector with filtered spikes - based on indexing from c
SpikeSpeed_fil=ST(c);
tSpk = SpikeSpeed_fil; % spike times

%% other stuff

% pull data out of struct
s = model.saved;
err = mean(s.FV, 2);
iter = s.IterCnt;
HD = get_hd(P);

% scale error (from 0-1)
err_flip = 1./err;
err_scaled = (err_flip-min(err_flip))*(1-0)/(max(err_flip)-min(err_flip)) + 0;

% parse fit parameters
g = s.ParamList(:,1);
theta_P = s.ParamList(:,2);
xref = s.ParamList(:,3);
yref = s.ParamList(:,4);

% calculate the 5th percentile
per5 = prctile(err,50); 
per5_idx = find(err<per5);
xref5 = xref(per5_idx);
yref5 = yref(per5_idx);
err5 = err(per5_idx);

% plot
hold on;
pathPlot_quiver(position, tSpk, get_hd(position));
% model_loc = scatter(xref, yref, [45], err_scaled, 'filled');
% map = brighten(flipud(hot),.65);
% colormap(map);
% model_loc.MarkerEdgeAlpha = 0.3;
best_fit = scatter(model.bestParams.xref, model.bestParams.yref, [120], 'MarkerEdgeColor','k', 'MarkerFaceColor', 'r', 'LineWidth', .8);
best_fit.MarkerFaceAlpha = 0.4;
real_loc = scatter(ref_point(1,1),ref_point(1,2),[100], 'MarkerEdgeColor','k', 'MarkerFaceColor', 'b', 'LineWidth', .8);
real_loc.MarkerFaceAlpha = 0.4;
pbaspect([1 1 1])%square
% hold off;

end

