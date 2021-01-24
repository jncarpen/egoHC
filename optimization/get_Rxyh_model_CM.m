function R_xyh_model = get_Rxyh_model_CM(pFit,rF,rP,rCutOff,Nbins)
% INPUT
%   rF:         R_data(x,y,H) (10x10x10)
%   rP:         r(x,y) (10x10)
%   rCutOff:    min firing rate for a bin to be considered
%   Nbins:      number of bins of the discretization
%   pFit:       [g; thetaP; xref; yref]; initial conditions 
%               of the 4 parameters to fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pull out parameters
g = pFit(1);
thetaP = pFit(2);
xref = pFit(3);
yref = pFit(4);

fF = 0;
    angCount = 0;
    for i=1:Nbins
        for j=1:Nbins
            if (rP(i, j)>rCutOff)
                rT = squeeze(rF(i, j, :)); % 10x1
                iF = find(isfinite(rT));
                if (length(iF)>0)
                    a = 180*atan2(yref-j, xref-i)/pi - (-180 -360/(2*Nbins) + iF*360/Nbins);
                    cFac = cos(pi*(a-thetaP)/180);
                    
                    % take circular mean- ignore NaN
                   cBar = circ_mean(cFac(~isnan(cFac)));
                    
                    z = 1+g*(cFac - cBar);
                    z = z.*(z>0);
                    R_xyh_model(i,j,:) = z;
                    fF = fF + nansum((z - rT(iF)).^2);
                    angCount = angCount + length(iF);
                end
            end
        end
    end
    fF = fF/angCount;
end