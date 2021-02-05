function [RXYHM fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins)
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

RXYHM = zeros(Nbins, Nbins, Nbins)*nan;

fF = 0;
    angCount = 0;
    for i=1:Nbins
        for j=1:Nbins
            if (rP(i, j)>rCutOff)
                rT = squeeze(rF(i, j, :));
                iF = find(isfinite(rT));
                iNF = find(~isfinite(rT));
                if ~isempty(iF)          
                    a = 180*atan2(yref-i, xref-j)/pi - ...
                        (-180 -360/(2*Nbins) + iF*360/Nbins);
                    cFac = cos(pi*(a-thetaP)/180);
                    cBar = nanmean(cFac);
                    
                    % shape the cosine
                    z = 1+g*(cFac - cBar);
                    z = z.*(z>0);
                    
                    % variance explained
                    fF = fF + nansum((z - rT(iF)).^2);                   
                    angCount = angCount + length(iF);
                    
                    % predicted output
                    RXYHM(i,j,iF) = z;
                    RXYHM(i,j, iNF) = nan;
                else
                    % predicted output
                    RXYHM(i,j,:) = nan;
                end
            end
        end
    end
    fF = fF/angCount;
end