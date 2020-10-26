% make box centers array

boxCtr = cell(1, length(pos_cm)); % initalize new cell array

for sessNum = 1:length(pos_cm)
    [boxCtrX,boxCtrY] = getBoxCenter(pos_cm{1,sessNum});
    
    boxCtr{1,sessNum}(1,1:2) = boxCtrX,boxCtrY;
end


% find home well locations for each session
hwCoord = cell(1,length(pos_cm));
for sessNum = 1:length(pos_cm)
    print_me = strcat(sessType{1,sessNum}, ' session #', sprintf('%.f', sessNum));
    disp(print_me)
    figure
    P = pos_cm{1,sessNum};
    plot(P(:,2), P(:,3))
    [x,y] = ginput(1);
    hwCoord{1,sessNum} = [x,y];
    
    close all
end
