function [corrPos] = Position_cm(pos, sessInfo, sessNum)
%POSITION_CM Summary of this function goes here
%   Correct position for a single session. This means convert it from bins
%   to centimeters and shift over so things should range around 0 to 150
%   cm since the box size is 1.5x1.5 m

% convert to centimeters and shift to start at 0
conversionFactor = 2.6;
minX = str2double(sessInfo{1,sessNum}.window_min_x{1,1});
minY = str2double(sessInfo{1,sessNum}.window_min_y{1,1});
posNow = pos{1,sessNum};

pos_ = zeros(length(posNow),5);
pos_(:,1) = posNow(:,1);
pos_(:,2) = (posNow(:,2)-minX)/conversionFactor;
pos_(:,3) = (posNow(:,3)-minY)/conversionFactor;
pos_(:,4) = (posNow(:,4)-minX)/conversionFactor;
pos_(:,5) = (posNow(:,5)-minY)/conversionFactor;

% remove outliers
x = pos_(:,2); y = pos_(:,3);
x2 = pos_(:,4); y2 = pos_(:,5);

d = [];
x_prime = []; y_prime = [];
x2_prime = []; y2_prime = [];

for idx = 1:length(x)-1
    d(idx) = sqrt((x(idx+1) - x(idx))^2 + (y(idx+1) - y(idx))^2);
    if d(idx) < 2
        x_prime(idx:idx+1) = x(idx:idx+1); y_prime(idx:idx+1) = y(idx:idx+1);
        x2_prime(idx:idx+1) = x2(idx:idx+1); y2_prime(idx:idx+1) = y2(idx:idx+1);
    else
        x_prime(idx:idx+1) = nan; y_prime(idx:idx+1) = nan;
        x2_prime(idx:idx+1) = nan; y2_prime(idx:idx+1) = nan;
    end
end

corrPos = [];
corrPos(:,1)=pos_(:,1); 
corrPos(:,2)=x_prime + abs(nanmin(x_prime)); corrPos(:,3)=y_prime + abs(nanmin(y_prime));
corrPos(:,4)=x2_prime + abs(nanmin(x2_prime)); corrPos(:,5)=y2_prime + abs(nanmin(y2_prime));

end

