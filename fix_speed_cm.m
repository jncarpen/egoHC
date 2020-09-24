function [speed_cm, accel_cm] = fix_speed_cm(pos_cm)
%FIX_SPEED_CM_FUNCT Summary of this function goes here
%   add a speed variable (speed in cm)
    totalSessions = length(pos_cm);
    speed_cm = cell(1,totalSessions);
    accel_cm = cell(1,totalSessions);

    for sessionNum = 1:totalSessions
        smoothPos = pos_cm{1,sessionNum};

        numLeds = 2; 
        % parse and medfilt pos vector to eliminate major outliers
        t = smoothPos(:,1); % in seconds
        x = medfilt1(smoothPos(:,2)); 
        x2 = medfilt1(smoothPos(:,4));
        y = medfilt1(smoothPos(:,3)); 
        y2 = medfilt1(smoothPos(:,5));
        v = zeros(size(smoothPos,1), numLeds); % velocity

        for i = 2:size(smoothPos,1)-1
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

        speed_cm{1,sessionNum} = v; 
        accel_cm{1,sessionNum} = a;
    end
    
end

