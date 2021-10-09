function [map] = raster_scan_integration(dx_map,dy_map,init_val)
%RASTER_SCAN_INTEGRATION Summary of this function goes here
%   Detailed explanation goes here
H = size(dx_map,1);
W = size(dy_map,2);
map = zeros(H,W);

map(1,1) = init_val;
for j=2:W
    map(1,j) = map(1,j-1) + dx_map(1,j-1);
end

for i=2:H
    for j=2:W
        map(i,j) = map(i-1,j) + dy_map(i-1,j);
    end
end
end

