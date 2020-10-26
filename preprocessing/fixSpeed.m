% Apply median filter to position and then re-calculate speed for all four
% LEDs to remove major outliers from the original speed vector. This fix
% has already been made in generateDataset.m so should not be a problem for
% the next animal...
%
% Jordan Carpenter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = pos_cm;

% initialize speed cell array
speed_cm = cell(1,length(P));

for sess = 1:length(P)
    
    % pull out current session
    t = P{1,sess}(:,1);
    x = medfilt1(P{1,sess}(:,2));
    y = medfilt1(P{1,sess}(:,3));
    x2 = medfilt1(P{1,sess}(:,4));
    y2 = medfilt1(P{1,sess}(:,5));
    v = zeros(length(t), 2);

    
    for i = 2:numel(x)-1
        v(i, 1) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
        v(i, 2) = sqrt((x2(i+1) - x2(i-1))^2 + (y2(i+1) - y2(i-1))^2) / (t(i+1) - t(i-1));
    end

    % pad the vector
    v(1,1) = v(2,1);
    v(1,2) = v(2,2);
    v(end, 2) = v(end-1, 2);
    v(end, 1) = v(end-1, 1);
    
    % save session to speed vector
    speed_cm{1,sess} = v;
       
end