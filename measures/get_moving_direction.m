function [md_1, md_2, md_3] = get_moving_direction(position)
%GET_MOVING_DIRECTION_TORGEIR Summary of this function goes here
%   Caclulate moving direction using different methods
%
%   Inputs:
%   'position'          position data (tx5 vector) for a single trial with this format:
%                       position = [t x1 y1 x2 y2] or position = [t x1 y1];
%
%   Outputs:            
%   'md_1'              tx1 vector, where t is the number of position
%                       samples, computed by using a space + time
%                       threshold.
%   'md_2'              same as md_1, but computing using the method
%                       Torgeir sent.
%   'md_3'              tx1 vector- based on the method Ulanovsky used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse and smooth position data
sigma = 3;
t = position(:,1); fs = mode(diff(t));
x = imgaussfilt(position(:,2), sigma);
y = imgaussfilt(position(:,3), sigma);
xy = [x,y];
md = diff(xy); % diff between consecutive tracking samples


%% METHOD 1

md_1 = [];
radius = 10; % cm
ffwd = floor(5/fs); % number of bins needed to ffwd 5 seconds
for tt = 1:length(t)
    if tt < length(t)-ffwd
        x_now = x(tt); y_now = y(tt);
        x_later = x(tt+ffwd); y_later = y(tt+ffwd); % position in 5 seconds
        d(tt) = sqrt((x_later-y_now).^2 + (y_later-y_now).^2); % distance between points
        
        if d(tt) >= radius
            md_1(tt) = rem(atan2d(md(tt,1), md(tt,2)) + 180, 360);
        else
            md_1(tt) = nan;
        end
        
    else
        % pad the vector with nans (for now)
        md_1(tt) = nan;
    end
end


%% METHOD 2
% Method Torgeir used:

% distance moved (hypotenuse between x y of samples)
movement_distance = sqrt(sum(md.^2, 2));

% direction moved (NaN when no movement)
md_2 = rem(atan2d(md(:,1), md(:,2)) + 180, 360);
md_2(movement_distance == 0) = nan;


%% METHOD 3
%   This method computes azimuth heading-direction (phi) based on the
%   following equation:
%
%               phi = angle(delta_x + delta_y * i)
%
%   where x and y are the position of the animal's head in space, delta_x
%   and delta_y are the changes in head-position between consecutive video
%   frames and the imaginary number i = sqrt(-1). This is the method used
%   in Sarel et al. 2017.
%   Something looks wrong here, need to look more into it (Oct 16, 2020)

clear i
md_3 = angle(md(:,1) + md(:,2) * 1i); % range is from -pi to +pi

% need to shift values

% convert to degrees (if you want)
md_3 = rad2deg(md_3+pi);
end

