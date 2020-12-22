% make a small population of simulated cells

%% egocentric place cell
sess=39;
param.position = pos_cm{1,sess};
param.theta = randi([0 360],1,1);
param.r = randi([10 40],1,1);

xRand = randi([ceil(nanmin(param.position(:,2))) floor(nanmax(param.position(:,2)))],1,2);
yRand = randi([ceil(nanmin(param.position(:,3))) floor(nanmax(param.position(:,3)))],1,2);

param.ctr_mass = [xRand(1), yRand(1)];
param.ref_point = [xRand(2) yRand(2)];

[sim] = simulate_place_ego(param);



%% egocentric bearing cell
clear param
sess = 50;
param.position = pos_cm{1,sess};
param.theta = randi([0 360],1,1);

xRand = randi([ceil(nanmin(param.position(:,2))) floor(nanmax(param.position(:,2)))],1,2);
yRand = randi([ceil(nanmin(param.position(:,3))) floor(nanmax(param.position(:,3)))],1,2);
param.ref_point = [xRand(2) yRand(2)];
param.r = randi([5 40],1,1);

[sim] = simulate_OVC(param);

%% hd cell
clear param
sess = 32;
param.position = pos_cm{1,sess};
param.theta = randi([0 360],1,1);
[sim] = simulate_HD(param);



%% assign
j = 12;
simdata(j).type = 'OVC';
simdata(j).spiketimes = sim.spiketimes;
simdata(j).position = param.position;
simdata(j).ctr_mass = [];%param.ctr_mass;
simdata(j).ref_point = param.ref_point;
simdata(j).theta = param.theta;

