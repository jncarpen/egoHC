function [MD] = get_MD(P)
%GET_MD 

t = P(:,1);
x = P(:,2);
y = P(:,3);

for ts = 1:length(t)-1  
    if ts == 1
        MD(1) = NaN;
    elseif ts == length(t)-1
        MD(end-1:end+2) = NaN;
    else
        % -pi to pi
        MD(ts) = atan2( y(ts+1) - y(ts-1), x(ts+1) - x(ts-1));
    end
end
end

