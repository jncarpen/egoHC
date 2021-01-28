function [R_xyh_model fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins)
% INPUT
%   rF:         R_data(x,y,H) (10x10x10)
%   rP:         r(x,y) (10x10)
%   rCutOff:    min firing rate for a bin to be considered
%   Nbins:      number of bins of the discretization
%   pFit:       [g; thetaP; xref; yref]; initial conditions 
%               of the 4 parameters to fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters to be fit
g = pFit(1);
thetaP = pFit(2); % (deg)
xref = pFit(3);
yref = pFit(4);

fF = 0;
    angCount = 0;
    for i=1:Nbins
        for j=1:Nbins
            if (rP(i, j)>rCutOff)
                rT = squeeze(rF(i, j, :));
                iF = find(isfinite(rT));
                iNF = find(~isfinite(rT));
                
                if (length(iF)>0)            
                    a = 180*atan2(yref-i, xref-j)/pi - ...
                        (-180 -360/(2*Nbins) + iF*360/Nbins);
                    
                   

                    cFac = cos(pi*(a-thetaP)/180);
                    cBar = nanmean(cFac);
                    
                    z = 1+g*(cFac - cBar);
                    % variance/error?
                    z = z.*(z>0);
                    fF = fF + nansum((z - rT(iF)).^2);                   
                    angCount = angCount + length(iF);
                    
%                     a_show = 180*atan2(yref-j, xref-i)/pi - ...
%                         (-180 -360/(2*Nbins) + [1:10]*360/Nbins);
%                     cFac_show = cos(pi*(a_show-thetaP)/180);
%                     cBar_show = nanmean(cFac_show);
%                     z_show = 1+g*(cFac_show - cBar_show);
                    R_xyh_model(i,j,iF) = z;
%                     R_xyh_model(i,j,iNF) = nan;
                end
            end
        end
    end
    fF = fF/angCount;
end

% pull out parameters
% g = pFit(1);
% thetaP = pFit(2);
% xref = pFit(3);
% yref = pFit(4);
% 
% R_xyh_model = zeros(10,10,10)*NaN;
% for i=1:Nbins
%     for j=1:Nbins
%         rT = squeeze(rF(i, j, :)); % 10x1
%         iF = linspace(1,10,10);
%         a = 180*atan2(yref-j, xref-i)/pi - (-180 -360/(2*Nbins) + iF*360/Nbins);
%         cFac = cos(pi*(a-thetaP)/180);
% %         cBar = nanmean(cFac);
%         cBar = circ_mean(cFac(~isnan(cFac)));
%         z = 1+g*(cFac - cBar);
%         z = z.*(z>0);
%         R_xyh_model(i,j,:) = z;
%     end
% end
    
% end