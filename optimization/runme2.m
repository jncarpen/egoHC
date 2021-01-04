%% load some sample data
load simdata

%% choose a cell
neuron = 4;
disp(strcat(simdata(neuron).type, ' cell'))

%% run the model
[model] = modelMe(simdata(neuron).position, simdata(neuron).spiketimes);

subplot(2,2,1)
plot_modelDynamics(simdata(neuron).position, simdata(neuron).spiketimes, model, simdata(neuron).ref_point)

subplot(2,2,2)
plot_vectorMod(model)

subplot(2,2,3)
plot_iterations(model)
