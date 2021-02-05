function [P_smooth] = smooth_pos(P_raw, sigma)
%SMOOTH_POS 
[~, c]=size(P_raw);
P_smooth(:,1) = P_raw(:,1); % we don't want to smooth time vector
for col = 2:c
    P_smooth(:,col) = imgaussfilt(P_raw(:,col), sigma);
end

end

