function [sim_pipj] = vsim(xi, xj, yi, yj, ui, uj, vi, vj)
%VECSIM: calculate similarity between two vectors

%% 1. alpha: euclidean distance
% calculate [Euclidean] distance term
dist_term = exp(-sqrt(((xi - xj).^2)+((yi - yj).^2)));

%% 2. beta: angular difference
veci = [ui vi];
vecj = [uj vj];

% take dot product of vi.vj
dot_vivj = dot(veci,vecj);

% calculate the norm of vi and vj
norm_vi = norm(veci);
norm_vj = norm(vecj);

% calcuate angular term
ang_term = exp(1-(dot_vivj./(norm_vi * norm_vj)));

%% 3. gamma: magnitude difference
% calculate magnitude term
mag_term = exp(-(norm(norm_vi-norm_vj)));

%% 4. vector pair similarity
% impose req. that all components contribute equally
alpha = 1/3;
beta = 1/3;
gamma = 1/3;

% calculate similarity value (range 0 to 1)
sim_pipj = alpha*dist_term + beta*ang_term + gamma*mag_term;
% sim_pipj = beta*ang_term + gamma*mag_term;


end

