function [valuesMapped] = map0to1(values)
%MAP0TO1 @jcarp2020

% check the input
% assert(isvector(vector), 'Please input a one-dimensional vector.')

% reshape the values
if size(values,1) > 1 && size(values,2) > 1
    vector = reshape(values, size(values,1)*size(values,1),1);
else
    vector = values;
end

% find min and max values of the input vector
min_vec = nanmin(vector);
max_vec = nanmax(vector);
span = max_vec-min_vec;

% convert vector into a 0-1 range
valuesMapped = (values - min_vec) ./ (span);

end
