function occupancyMap(posx, posy, speed_s)
%OCCUPANCYMAP: create occupancy and speed maps

% plot occ and map
for i=1:100
for j= 1:100
f = find(abs(posx-i)<5 & abs(posy-j)<5);
occ(i,j) = length(f);
map(i,j) = nanmedian(speed_s(f));
end
end

figure 
imagesc(occ)
title("Occupancy")
xlabel("X-Coordinate")
ylabel("Y-Coordinate")
colorbar

figure
contourf(occ)
title("Occupancy- Contour")
xlabel("X-Coordinate")
ylabel("Y-Coordinate")
colorbar

figure
imagesc(map)
title("Speed Map")
xlabel("X-Coordinate")
ylabel("Y-Coordinate")
colorbar

figure
contourf(map)
title("Speed Map- Contour")
xlabel("X-Coordinate")
ylabel("Y-Coordinate")
colorbar

return
end

