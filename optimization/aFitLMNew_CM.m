function [fF] = aFitLMNew_CM(pFit,rF,rP,rCutOff,Nbins)
%   INPUT -
%   rF:             R_data(x,y,H) (10x10x10)
%   rP:             r(x,y) (10x10)
%   rCutOff:        min firing rate for a bin to be considered
%   Nbins:          number of bins of the discretization
%   pFit:           [g; thetaP; xref; yref]; initial conditions 
%                   of the 4 parameters to fit
%   OUTPUT -
%   fF:             error (?)
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
                % angular tc for bin(i,j)
                rT = squeeze(rF(i, j, :));
                
                % index of finite values in tc (ignores NaN & inf)
                iF = find(isfinite(rT));
                
                if (length(iF)>0)
                    
                   % a = 180*atan2(yref-j, xref-i)/pi - (-180 -360/(2*Nbins) + iF*360/Nbins);
                   
                   % find the allocentric bearing (-180:180 deg)
                   allo = atan2d(yref-j, xref-i);
                   
                   % calculate shift
                   shift = -180 -360/(2*Nbins) + iF*360/Nbins;
                   
                   % shift (this calculates ego bearing)
                   theta = mod(allo - shift, 360)-180;
                   
                   % cosine term 
                   cFac = cos(deg2rad(theta-thetaP));
                   
                   % take the mean of the cosine (circular)
                   cBar = circ_mean(cFac(~isnan(cFac))); 
                   
                   % shift and stretchs
                   z = 1+g*(cFac - cBar);
                   
                   % only keep nonzero elements
                   z = z.*(z>0);
                   
                   % error
                   fF = fF + nansum((z - rT(iF)).^2);
                   angCount = angCount + length(iF);
                    
                end
            end
        end
    end
    fF = fF/angCount;
end