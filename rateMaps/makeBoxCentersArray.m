boxCtr = cell(1, length(pos_cm)); % initalize new cell array

for sessNum = 1:length(pos_cm)
    [boxCtrX,boxCtrY] = getBoxCenter(pos_cm{1,sessNum});
    
    boxCtr{1,sessNum}(1,1:2) = boxCtrX,boxCtrY;
end
