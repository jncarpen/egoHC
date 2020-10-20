boxSize = 80; % for Sebastian's data
[pos_in] = correct_pos_general(positions, boxSize);
ref_point = [42, 38];
AOI = 15:15:345;
ROI = 15:5:40;

for aa = 1:length(AOI)
    for rr = 1:length(ROI)
        % make a new figure
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color','w');
        
        % get values for this iteration
        angle_of_interest = AOI(aa);
        radius = ROI(rr);

        % set figure title
        figTit = strcat('preferred angle: ', sprintf('%.f', angle_of_interest), 'deg, radius:', sprintf('%.f', radius), 'cm');
        sgtitle(figTit);


        % egocentric bearing
        subplot(2, 3, 1)
        [~, ~, ~] = simulate_ego_cell(pos_in, ref_point, angle_of_interest);

        % egocentric bearing + distance
        subplot(2, 3, 2)
        [~, ~, ~] = simulate_egoDist_cell(pos_in, ref_point, angle_of_interest, radius);

        % allocentric bearing
        subplot(2, 3, 3)
        [~, ~, ~] = simulate_allo_cell(pos_in, ref_point, angle_of_interest);

        % allocentric bearing + distance (OVC)
        subplot(2, 3, 4)
        [~, ~, ~] = simulate_OVC(pos_in, ref_point, angle_of_interest, radius);

        % head-direction cell
        subplot(2, 3, 5)
        [~, ~, ~] = simulate_HD(pos_in, angle_of_interest);

        % save figure
        filename = strcat('D:\egoAnalysis\simulated_cells_coverage\', 'PA ', sprintf('%.f', angle_of_interest), 'R', sprintf('%.f', radius), '.png');
        saveas(fig, filename);
        
        close all; clf;
        
    end
end
