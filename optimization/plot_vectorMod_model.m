function plot_vectorMod_model(model)
clear bins_corr

% number of bins used
nBins = 10;
% calculate binwise error. @todo: make small values map onto big ones
for row = 1:nBins
    for col = 1:nBins
        if ~isempty(model.predcell{row,col}) && ~isempty(model.datacell{row,col})
            binwise_err(row,col) = sum((model.predcell{row,col}(:)-model.datacell{row,col}(:)).^2);
            % find peak orientations of tuning curves
            [~, ipred] = nanmax(model.predcell{row,col}(:));
            peak_pred(row,col) = model.bins(ipred);
            [~, idata] = nanmax(model.datacell{row,col}(:));
            peak_data(row,col) = model.bins(idata);
            % find the mean vector length of data
            data_tcStats = analyses.tcStatistics(model.datacell{row,col}(:), mode(diff(model.bins)), 95);
            data_MVL(row,col) = data_tcStats.r; 
        else
            binwise_err(row,col) = NaN; peak_pred(row,col)= NaN; peak_data(row,col)=NaN; data_MVL(row,col)=NaN;
        end
    end
end

% reshape vectors, @todo use phase angle instead of peak?
pred_val = reshape(peak_pred', 100, 1);
data_val = reshape(peak_data', 100, 1);

% reshape MVL (scaling factor)
data_MVL = reshape(model.modStrength.HD_MVL, 100, 1);
data_MVL(isnan(data_MVL))=0;

model_MVL = reshape(model.modStrength.RH_MVL, 100, 1);
model_MVL(isnan(model_MVL))=0;

% map small values onto big values
% fac = 2;
% fac_data =.75;
% pred_scale_unmapped = reshape((1./binwise_err)', 100, 1);
% % map scale for predicted data onto 0-1
% pred_scale = (pred_scale_unmapped-min(pred_scale_unmapped))*(1-0)/(max(pred_scale_unmapped)-min(pred_scale_unmapped)) + 0;
% pred_scale = pred_scale.*fac;
% data_scale = reshape(data_MVL', 100, 1).*fac_data;

% define vector orientations
u = cos(pred_val * pi/180); 
v = sin(pred_val * pi/180);
u_data = cos(data_val * pi/180);
v_data = sin(data_val * pi/180);

% find scaling factor
% sf = abs(pred_scale./(sqrt((u.^2)+(v.^2))));
% sf_data = abs(data_scale./(sqrt((u_data.^2)+(v_data.^2))));

sf = model_MVL;
sf_data = data_MVL;

% scale the vectors
uprime = u.*sf; 
vprime = v.*sf;
uprime_data = u_data.*sf_data; 
vprime_data = v_data.*sf_data;

% get rid of nan values in ratemap
rm_vec = reshape(model.rateMap',100,1);
nan_idx = find(isnan(rm_vec));
bins_corr = model.spatbinsnum; % make copy
bins_corr(nan_idx,:) = NaN;

% figure
set(gca, 'visible', 'off')
hold on;
% imagescwithnan(flipud(model.rateMap),jet,[1 1 1])
% check if 
% imagescwithnan(flipud(model.rateMap),jet,[.7 .5 .7])% colorbar
% alpha(0.3) 
% brighten(.6)
xlim([0 11]); ylim([0 11]);
% plot data vectors
% modelVecs = quiver(bins_corr(:,1), bins_corr(:,2), uprime_data, vprime_data, 0);
% set(modelVecs, 'Color', 'blue', 'AutoScale', 'off', 'LineWidth',1)
% plot model vectors
modelVecs2 = quiver(bins_corr(:,1), bins_corr(:,2), uprime, vprime, 0);
set(modelVecs2, 'Color', 'r', 'AutoScale', 'off', 'LineWidth',1)
% make the figure a square
pbaspect([1 1 1])

end




% x_scaled = (x-min(x))*(b-a)/(max(x)-min(x)) + a;
% pred_val = repmat(model.thetaP, 100,1);
% apply an offset value (unnecessary)
% offset = atan2d(model.yref-model.spatbins(:,2), model.xref-model.spatbins(:,1))+180;
