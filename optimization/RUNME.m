%% December 7, 2020 RUN ME

%% pick a session number
sessNum = 1;

%% (1) regular place cell
clear param sim
% set parameters for simulated cell
param.position = data(sessNum).pos;
param.ctr_mass = [75 75]; % center of arena
param.noise = 0; % no noise 
param.width = 5;
% simulate the cell
[sim] = simulate_place(param);
% run the model
[model] = modelMe(sim.position, sim.spiketimes);
ref_point = [model.bestParams.xref, model.bestParams.yref];
% visualize
subplot(1,3,1)
plot_modelDynamics(sim.position, sim.spiketimes, model, ref_point);
subplot(1,3,2)
plot_vectorMod(model)
subplot(1,3,3)
plot_iterations(model)



%% directionally-modulated place cell
clear param sim
% set parameters for simulated cell
param.position = data(sessNum).pos;
param.ctr_mass = [75 75]; % center of arena
param.noise = 0; % no noise 
param.width = 5;
param.theta = 270;
% simulate the cell
[sim] = simulate_place_egoMod(param);
% run the model
[model] = modelMe(sim.position, sim.spiketimes);
ref_point = [model.bestParams.xref, model.bestParams.yref];
% visualize
subplot(1,3,1)
plot_modelDynamics(sim.position, sim.spiketimes, model, ref_point);
subplot(1,3,2)
plot_vectorMod(model)
subplot(1,3,3)
plot_iterations(model)

