% Test to see how the BNT script works

addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001\bnt"));

for sess = 12
    P = pos{1,sess};
    PXY = pos{1,sess}(:,2:3);
    ST = SpikeTimes{1,sess};
    for unit = 1:length(ST)
        STC = ST{1,unit};
        % map = analyses.map(P, STC, 'BinWidth',4);
        [spkPos, spkInd] = data.getSpikePositions(STC, P);
        [map, posPdf] = analyses.mapAdaptiveSmoothing(PXY, spkPos);
        % peak = nanmax(nanmax(map.z));
        
        figure
        % subplot(2,2,1)
        imagesc(map.z)
        title("z")
        colorbar
        % caxis([.2 peak])
        
%         subplot(2,2,2)
%         imagesc(map.count)
%         title("SpikeCount")
%         colorbar
%         
%         subplot(2,2,3)
%         imagesc(map.time)
%         title("occupancy")
%         colorbar
    end
end