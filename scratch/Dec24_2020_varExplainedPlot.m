%% December 24, 2020
% Check to see if all variance explained values are between 0 and 1 and if 
% they can be plotted as they are in the Jercog paper

clear varex_place varex_model modstren_hd modstren_rh err
count = 1;
for nn = 1000:1500 %length(JZ.neurons)
    for uu = 1:length(JZ.neurons(nn).members)
        
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(nn).members(uu).P);
        ST = JZ.neurons(nn).members(uu).ST;

        % smooth position vectors
        clear P
        sigma = 2;
        P = smooth_pos(P_raw, sigma);
        HD = get_hd(P);
        
        initial = choose_initial_conditions(P);
        [model] = modelMe(P, ST, HD, initial);
        
        varex_place(count) = model.varExplained.place;
        varex_model(count) = model.varExplained.model;
        
        modstren_hd(count) = model.modStrength.HD;
        modstren_rh(count) = model.modStrength.RH;
        
        err(count) = model.err;
        
        
        count = count + 1;

    end
end



% name the variables
X = varex_model;
plotName = 'variance explained by model';
xname = 'variance explained';

% plot
figure; set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceColor', [1 .2 0]);
ylabel("frequency (count)"); xlabel(xname);
title(plotName);
set(gca,'FontSize',20, 'FontName', 'Calibri Light', 'FontWeight', 'normal');
box off;

