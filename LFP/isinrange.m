function [v,ix,iw,n] = isinrange(x, xw, mode, check)
%ISINRANGE returns the indices of the sorted numeric vector x that fall within the
%windows with initial positions designated by xw and lengths specified by
%the wsz. N.B. this will not work if X is unsorted; use isinrange2 instead!
% 
% INPUTS
% x    - data vector
% xw   - a nwin*2 matrix of 'start' and 'stop' values for each window.
% wsz  - in the case where xw is a vector, a scalar indicating window size,
%
% OUTPUTS
% v     - logical vector of with length equal to length(x).  For x(i), v(i)
%         is true if x(i) occurs within a window, and is false otherwise.
% ix    - nwindows*2 matrix of index ranges of x in each window, such that
%         x(ix(7,1) : ix(7,2)) returns all x values within the 7th window
% iw    - vector of window indices for each value in x
% n    -  1*nwindows vector of number of values of x within each window
% 
% 
% R Gardner 2013-03-29
% 2013-06-22 add output for n
% 2019-01-31 check for sorted input by default
% 2019-03-26 

if nargin < 4 || isempty(check), check = true; end

v = false(size(x));
iw = nan(size(x));

if iscolumn(x)
    x = x';
end

nans = isnan(x);
xn = x(~nans);
xni = find(~nans);

if check
    assert( ...
        issorted(xn), ...
        'rg:isinrange:xNotSorted', ...
        'Input "x" must be a sorted numeric vector.');
end

if size(xw,2) ~= 2
    error('xw must be a 2-column matrix')
end

if nargin < 3 || isempty(mode)
    mode = 'incl';
end

nw = size(xw,1);
nxn = numel(xn);

if isempty(xw) || isempty(x)
    ix = [];
    n = zeros(nw, 1);
    return
end

% Find the index for the first and last x value before each win start/stop
i1 = rg.helpers.m_lookup(xw(:,1)', xn); % Last BEF / EQUAL TO start
i2 = rg.helpers.m_lookup(xw(:,2)', xn); % Last BEF / EQUAL TO end

% if xw(end, 1), its index will be zero; correct this
% to nx
i1(xw(:, 1)==xn(end)) = nxn;

% Correct not-found indices
i2(i2==0 & i1>0) = nxn;  % Only fix end indices when a start index WAS found

% Special case: single window where all values of x lie within bounds
if nw==1 && i2==0 && i1==0
    if xw(1) <= xn(1) && xw(2) >= xn(end)
        i2 = nxn;
    end
end

% Change i1 so that it contains the index to first x value AFTER win start
i1 = i1+1;

if strcmpi(mode,'incl')
    
    % Find x values that match window starts and shift their indices one to
    % the left
    match1 = ismembc(xw(:,1), xn);
    i1(match1) = i1(match1)-1;
    
    if xn(1)==xw(1,1)
        i1(1) = 1;
    end
    
elseif strcmpi(mode,'excl')
    % Find x values that match window ends and shift their indices one to
    % the left    
    match2 = ismembc(xw(:,2), xn);
    i2(match2) = i2(match2)-1;
else
    error('Invalid inclusivity mode')
end

% Convert x back to original form including NaN values
v1 = i1 ~= 0;
i1(v1) = xni(i1(v1));
v2 = i2 ~= 0;
i2(v2) = xni(i2(v2));

% Loop through windows
for i=1:nw
    % update logical output
    v(i1(i):i2(i)) = true;
    iw(i1(i):i2(i)) = i;
end
    
if nargout > 1
    n = i2'-i1'+1;
    ix = [i1' i2'];
end

end

