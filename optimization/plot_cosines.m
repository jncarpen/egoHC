function plot_cosines(model)


for row = 1:10
    for col = 1:10
        plot(model.bins, model.predcell{row,col})
        hold on
    end
end


end

