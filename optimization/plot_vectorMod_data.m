function plot_vectorMod_data(model)
% reshape vectors, @todo use phase angle instead of peak?
% +180 brings 'forward' to 'up' (not sure what to do here)
data_val = reshape(model.modStrength.HD_prefVec, 100, 1);

% reshape MVL (scaling factor)
data_MVL = reshape(model.modStrength.HD_MVL, 100, 1);
data_MVL(isnan(data_MVL))=0;

% define vector orientations
u_data = cos(data_val * pi/180);
v_data = sin(data_val * pi/180);

% scaling factors
sf_data = data_MVL;

% scale the vectors
uprime_data = u_data.*sf_data; 
vprime_data = v_data.*sf_data;

% get rid of nan values in ratemap
rm_vec = reshape(model.rateMap',100,1);
nan_idx = find(isnan(rm_vec));
binX = reshape(model.spatBinNum.X, 100, 1);
binY = reshape(model.spatBinNum.Y, 100, 1);


%% PLOT
figure; hold on;
set(gca, 'visible', 'off')

% plot data vectors
modelVecs = quiver(binX, binY, uprime_data, vprime_data, 0);
set(modelVecs, 'Color', 'blue', 'AutoScale', 'off', 'LineWidth',1)

% make the figure a square
pbaspect([1 1 1])

end


