%% December 22, 2020
% Define the false positive rate of the model
% Using the JZ.mat data
% J. Carpenter

%% (1) SHUFFLE HD 
% Create surrogate data in which the HD for the animal was shifted in 
% time relative to the location and firing rate. We shifted in a circular 
% manner with the end of the session wrapped to the beginning.


%% (2) SURROGATE HD
% Create surrogate data by randomly interchanging indices of the HD angles 
% at each time (draw without replacement). Repeat this procedure 1000 times
% and obtain a 95% CI limits and mean for the surrogated HD tuning curves.
