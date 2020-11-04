function err= WebersLawErr(p,x,y)

% DESCRIPTION: Here's a two-line function 'WebersLawErr' that calculates 
% the sums of squared error between the model's prediction and the data.

pred = WebersLaw(p,x); 
err = sum( (pred(:)-y(:)).^2);

end