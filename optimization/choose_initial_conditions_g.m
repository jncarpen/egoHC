function initial = choose_initial_conditions_g()
%CHOOSE_INITIAL_CONDITIONS Summary of this function goes here

% randomly sample from them 
howMany = 1;
initial.g = rand(howMany, 1);
% orientation = 1:1:360;
initial.thetaP = 0;%randsample(orientation,howMany)';
% initial.thetaP = 0;
initial.xref = 32;
initial.yref = 26;
end

