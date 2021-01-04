function initial = choose_initial_conditions(P)
%CHOOSE_INITIAL_CONDITIONS Summary of this function goes here

% make a vector of all possible positions
step = 0.5;
x = min(P(:,2)):step:max(P(:,2));
y = min(P(:,3)):step:max(P(:,3));
orientation = 1:1:360;

% randomly sample from them 
howMany = 1;
initial.g = rand(howMany, 1).*2.5;
initial.thetaP = randsample(orientation,howMany)';
initial.xref = randsample(x,howMany)';
initial.yref = randsample(y,howMany)';
end

