function plot_vectorMod_model(model)
% reshape vectors, @todo use phase angle instead of peak?
% +180 brings 'forward' to 'up' (not sure what to do here)
pred_val = reshape(model.modStrength.RH_prefVec, 100, 1);
data_val = reshape(model.modStrength.HD_prefVec, 100, 1);

% reshape MVL (scaling factor)
data_MVL = reshape(model.modStrength.HD_MVL, 100, 1);
data_MVL(isnan(data_MVL))=0;

model_MVL = reshape(model.modStrength.RH_MVL, 100, 1);
model_MVL(isnan(model_MVL))=0;

% define vector orientations
u = cos(pred_val * pi/180); 
v = sin(pred_val * pi/180);
u_data = cos(data_val * pi/180);
v_data = sin(data_val * pi/180);

% scaling factors
sf = model_MVL.*3;
sf_data = data_MVL;

% scale the vectors
uprime = u.*sf; 
vprime = v.*sf;
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
xlim([0 11]); ylim([0 11]);

% plot model vectors
modelVecs = quiver(binX, binY, uprime, vprime, 0);
set(modelVecs, 'Color', 'r', 'AutoScale', 'off', 'LineWidth',1)

% make the figure a square
pbaspect([1 1 1])
end




% x_scaled = (x-min(x))*(b-a)/(max(x)-min(x)) + a;
% pred_val = repmat(model.thetaP, 100,1);
% apply an offset value (unnecessary)
% offset = atan2d(model.yref-model.spatbins(:,2), model.xref-model.spatbins(:,1))+180;
