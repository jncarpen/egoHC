function [posCm, velocity_cm] = pix2cm(pos, sessInfo)
%PIX2CM Summary of this function goes here
%   Assuming the box is 1x1m in size.
%   1m = 100cm
%   I am now learning that I might not need this function at all...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% break up the pos matrix
t = pos(:,1);
x = medfilt1(pos(:,2));
y = medfilt1(pos(:,3));
x2 = medfilt1(pos(:,4));
y2 = medfilt1(pos(:,5));

% grab max/min x/y values from sessInfo struct
max_x = str2double(sessInfo.window_max_x{1,1});
min_x = str2double(sessInfo.window_min_x{1,1});
max_y = str2double(sessInfo.window_max_y{1,1});
min_y = str2double(sessInfo.window_min_y{1,1});

% compute window range for x/y
Xfactor = 100/(max_x - min_x);
Yfactor = 100/(max_y - min_y);

% create a new position matrix
posCm = [t, Xfactor.*x, Yfactor.*y, Xfactor.*x2, Yfactor.*y2];


%% Compute velocity given position in cm/s

numLEDs = 2; % can incorporate this into a diff function or make 
             % it toggle-able for this fn.

% redefine x/y          
t_cm = pos(:,1);
x_cm = pos(:,2);
y_cm = pos(:,3);
x2_cm = pos(:,4);
y2_cm = pos(:,5);

% initialize velocity vector
velocity_cm = zeros(size(posCm,1), numLEDs);

for i = 2:size(posCm,1)-1
    velocity_cm(i, 1) = sqrt((x_cm(i+1) - x_cm(i-1))^2 + (y_cm(i+1) - y_cm(i-1))^2) / (t_cm(i+1) - t_cm(i-1));
    velocity_cm(i, 2) = sqrt((x2_cm(i+1) - x2_cm(i-1))^2 + (y2_cm(i+1) - y2_cm(i-1))^2) / (t_cm(i+1) - t_cm(i-1));
end

% pad velocity vector
velocity_cm(1,1) = velocity_cm(2,1);
velocity_cm(1,2) = velocity_cm(2,2);
velocity_cm(end, 2) = velocity_cm(end-1, 2);
velocity_cm(end, 1) = velocity_cm(end-1, 1);

end

