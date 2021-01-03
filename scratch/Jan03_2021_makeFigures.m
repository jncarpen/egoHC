% make figures

% load dataset from the first animal
load("D:\Data\Dataset\24116\24116_2.mat")

% choose position information from a random session
P = pos_cm{1, randi(length(pos_cm))};
possible_x = min(P(:,2))+10:step:max(P(:,2))-10;
possible_y = min(P(:,3))+10:step:max(P(:,3))-10;

% set parameters for simulation
clear param
param.position = P;
param.ctr_mass = [rand(1)*150, rand(1)*150];
param.noise = 0;
param.width = rand(1)*10;

% simulate a place cell
[sim] = simulate_place(param);



