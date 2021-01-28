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
R_xyh_model = zeros(10,10,10)*NaN;
    for i=1:Nbins
        for j=1:Nbins
            rT = squeeze(rF(i, j, :)); % 10x1
            iF = linspace(1,10,10);
            a = 180*atan2(yref-j, xref-i)/pi - (-180 -360/(2*Nbins) + iF*360/Nbins);
            cFac = cos(pi*(a-thetaP)/180);
            cBar = circ_mean(cFac(~isnan(cFac)));
            z = 1+g*(cFac - cBar);
            z = z.*(z>0);
            R_xyh_model(i,j,:) = z;
        end
    end

end