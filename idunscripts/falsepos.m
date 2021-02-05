function falsepos(jobnum, filepath, datacodepath)
%FALSEPOS 

% add path where code is being stored
addpath(genpath(datacodepath));

% load data
load(strcat(datacodepath, '\', 'placeCell.mat'));

whichUnit = randi(100);
now = unimodal(whichUnit);
P = placeCell(now).P;
ST = placeCell(now).ST;
Z = placeCell(now).HD;

% make spiketrain
t = P(:,1); tpf = mode(diff(t));
startTime = t(1); stopTime = t(end);
ST = ST(ST < stopTime & ST > startTime);
t_edges = linspace(startTime,stopTime,numel(t)+1);
SpkTrn = histcounts(ST,t_edges);

% make ratemap
bins = 20;
[~, ~, ~, binX, binY] = histcounts2(P(:,2),P(:,3),bins);
binX(binX==0)=nan; binY(binY==0)=nan;
r_xy = zeros(bins,bins);
for i = 1:bins
    for j = 1:bins
        idx_here = find(i == binY & j == binX);
        spikes_here = SpkTrn(idx_here);
        if ~isempty(idx_here)
            r_xy(i,j) = sum(spikes_here)/(length(idx_here)*tpf);
        end
    end
end
map.z = flipud(r_xy);
map.whichBin.x = binX;
map.whichBin.y = binY;

% simulate place cell
[sim] = simulate_place(map, P(:,1));
ST = sim.ST;

% save some stuff
FP.ST = ST;
FP.param = root;

% sampling frequency (50 samples/sec)
Fs = mode(diff(P(:,1))); 

% perform the optimization (using a monte carlo method)
total_iters = 100; 
clear monte error_for_comparison
for optim_iter = 1:total_iters 
    % run the model 
    out = modelMe(P, ST, Z);
    % save all the runs to be compared
    monte(optim_iter).model = out; 
    error_for_comparison(optim_iter) = out.model.error;
end

% find run that yieled smallest error value
[~, errorValsMin] = nanmin(error_for_comparison);

% save the run that gives the global min
modelData = monte(errorValsMin).model;

% save stuff
FP.modelData = modelData;


%% shuffle the head direction values
% define the number of shuffles we want
total_shuffles = 1000;
min_shuffle = floor(60/Fs); % 30 seconds
max_shuffle = floor(600/Fs); % 200 seconds

% systematic shuffle (can be replaced with a pseudorandom shuffle)
shuffle_vec = [linspace(-min_shuffle, -max_shuffle, total_shuffles/2), ...
    linspace(min_shuffle, max_shuffle, total_shuffles/2)]; 
shuffle_vec =  floor(shuffle_vec); % make them all integers

% shuffle the head direction values 1000x
for shuff_num = 1:total_shuffles
    clear model_shuffled
    % grab the current shuffle value from vector
    shuff_value = shuffle_vec(shuff_num);

    % circularly shifts the elements in array A by K positions
    shuffled_Z = circshift(Z, shuff_value);

    % run the model on the shuffled data
    [out_SH] = modelMe(P, ST, shuffled_Z);

    % SAVE OUTPUT FROM EACH RUN

    % (1) best fit parameters
    param_g(shuff_num) = out_SH.model.fitParams.g;
    param_thetaP(shuff_num) = out_SH.model.fitParams.thetaP;
    param_xref(shuff_num) = out_SH.model.fitParams.xref;
    param_yref(shuff_num) = out_SH.model.fitParams.yref;

    % (2) variance explained (place & model)
    ve_place(shuff_num) = out_SH.measures.VE.place;
    ve_rh(shuff_num) = out_SH.measures.VE.RH;

    % (3) modulation/tuning strength
    TS_hd(shuff_num) = out_SH.measures.TS.HD;
    TS_rh(shuff_num) = out_SH.measures.TS.RH;

    % (4) error
    error_shuff(shuff_num) = out_SH.model.error;

end

% save stuff
FP.ve_place = ve_place;
FP.ve_rh = ve_rh;
FP.ts_hd = TS_hd;
FP.ts_rh = TS_rh;
FP.shuffErr = error_shuff;
FP.g = param_g;
FP.thetap = param_thetaP;
FP.xref = param_xref;
FP.yref = param_yref;

% save
fileBody = strcat('job', sprintf('%.f', jobnum));
savename = strcat(filepath, '\', fileBody, '.mat');
save(savename, 'FP', '-v7.3');


end

