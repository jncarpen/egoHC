%% December 22, 2020
% Define the false positive rate of optimization procedure
% Using the JZ.mat data (egoHC dataset from dMan)
% J. Carpenter


%% load data/ add path 
% load("D:\Data\Dataset\dManCopyDataV1\JZ.mat") % loads in JZ.mat file
% addpath(genpath("C:\Users\17145\Documents\github_local"))
% addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001\bnt"))


%% clear stuff
clear optimSig

%% run the model for one of the units
% choose a neuron/unit to check
% nn = 587; 
% uu = 1;
% 
% % pull information for this neuron
% P_raw = units2cm(JZ.neurons(nn).members(uu).P);
% ST = JZ.neurons(nn).members(uu).ST;
% Fs = mode(diff(P(:,1))); % sampling frequency (50 samples/sec)
% 
% % smooth position vectors
% sigma = 2; % width of Gaussian kernel
P = smooth_pos(P_raw, sigma);

% run the model on a simulated cell first
simParams.position = P; simParams.ref_point = [75, 75];
simParams.theta = 270;
[sim] = simulate_ego_cell(simParams);

% get 'real' head direction values (deg)
HD = get_hd(P);

% perform the optimization (using a monte carlo method)
total_iters = 100; 
clear monte error_for_comparison
for optim_iter = 1:total_iters 
    % choose some initial conditions randomly
    initial = choose_initial_conditions(P);
    % run the model 
    [model] = modelMe(P, ST, HD, initial);
    
    monte(optim_iter).model = model; % save all the runs to be compared
    error_for_comparison(optim_iter) = model.err;
end

% find run that yieled smallest error value
[~, errorValsMin] = nanmin(error_for_comparison);

% save the run that gives the global min
% optimSig(nn).members(uu).modelReal = monte(errorValsMin).model;
optimSig.modelReal = monte(errorValsMin).model;


%% shuffle the head direction values
% define the number of shuffles we want
total_shuffles = 100;
min_shuffle = floor(15/Fs); % 15 seconds
max_shuffle = floor(120/Fs); % 120 seconds

% systematic shuffle (can be replaced with a pseudorandom shuffle)
shuffle_vec = [linspace(-min_shuffle, -max_shuffle, total_shuffles/2), ...
    linspace(min_shuffle, max_shuffle, total_shuffles/2)]; 
shuffle_vec =  floor(shuffle_vec); % make them all integers

% clear stuff from before
clear param_g param_thetaP param_xref param_yref ...
    varex_place varex_rh modstren_hd modstren_hd error_shuff

% shuffle the head direction values 100x
for shuff_num = 1:total_shuffles
    clear model_shuffled
    % grab the current shuffle value from vector
    shuff_value = shuffle_vec(shuff_num);
    
    % circularly shifts the elements in array A by K positions
    shuffled_hd = circshift(HD, shuff_value);

    % generate random values for parameter intial conditions
    initial = choose_initial_conditions(P);

    % run the model on the shuffled data
    [model_shuffled] = modelMe(P, ST, HD, initial);
    
    % save data from each run
    % (1) best fit parameters
    param_g(shuff_num) = model_shuffled.bestParams.g;
    param_thetaP(shuff_num) = model_shuffled.bestParams.thetaP;
    param_xref(shuff_num) = model_shuffled.bestParams.xref;
    param_yref(shuff_num) = model_shuffled.bestParams.yref;
    % (2) variance explained (place & model)
    varex_place(shuff_num) = model_shuffled.varExplained.place;
    varex_rh(shuff_num) = model_shuffled.varExplained.model;
    % (3) modulation/tuning strength
    modstren_hd(shuff_num) = model_shuffled.modStrength.HD;
    modstren_rh(shuff_num) = model_shuffled.modStrength.RH;
    % (4) error
    error_shuff(shuff_num) = model_shuffled.err;
end

% name the variables
X = modstren_rh;
Xd = optimSig.modelReal.modStrength.RH;
plotName = 'modulation strength (rh)';
xname = 'modulation strength';

% plot
figure; set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', [0 .8 .8]);
xline(Xd, '--k', 'LineWidth', 1.5);
ylabel("frequency (count)"); xlabel(xname);
title(plotName);
set(gca,'FontSize',20, 'FontName', 'Calibri Light', 'FontWeight', 'normal');
l = legend('shuffled data', 'real data');
box off;

%% plot stuff about the current cell
pathPlot_hd(P, ST, HD)








