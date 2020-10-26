function [speed, accel] = get_speed(position)
%FIX_SPEED_CM_FUNCT Summary of this function goes here
%   get speed & acceleration

    numLeds = 2; 
    % parse and medfilt pos vector to eliminate major outliers
    t = position(:,1); % in seconds
    x = medfilt1(position(:,2)); 
    x2 = medfilt1(position(:,4));
    y = medfilt1(position(:,3)); 
    y2 = medfilt1(position(:,5));
    v = zeros(size(position,1), numLeds); % velocity

    for i = 2:size(position,1)-1
        v(i, 1) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
        v(i, 2) = sqrt((x2(i+1) - x2(i-1))^2 + (y2(i+1) - y2(i-1))^2) / (t(i+1) - t(i-1));
    end

    % pad the vector
    v(1,1) = v(2,1);
    v(1,2) = v(2,2);
    v(end, 2) = v(end-1, 2);
    v(end, 1) = v(end-1, 1);

    v1 = v(:,1);
    v2 = v(:,2);

    for i = 2:size(v,1)-1
      a(i,1) = (v1(i+1) - v1(i-1))/(t(i+1) - t(i-1));
      a(i,2) = (v2(i+1) - v2(i-1))/(t(i+1) - t(i-1));
    end

    speed = v;
    accel = a;
end



