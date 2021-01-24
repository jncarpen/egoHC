function initial = choose_initial_conditions(nbins)
%CHOOSE_INITIAL_CONDITIONS Summary of this function goes here

% make a vector of all possible positions
step = 0.5;
x = 1:5.:nbins;
y = 1:.5:nbins;
orientation = -180:1:180;

% randomly sample from them 
howMany = 1;
initial.g = rand(howMany, 1).*2.5;
initial.thetaP = randsample(orientation,howMany)';
initial.xref = randsample(x,howMany)';
initial.yref = randsample(y,howMany)';
end

