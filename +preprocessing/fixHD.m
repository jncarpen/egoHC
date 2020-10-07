newHD = cell(1,length(hd));

for i = 1:length(hd)
    curr = hd{1,i};
    newHD{1,i} = inpaint_nans(curr,3);
end
