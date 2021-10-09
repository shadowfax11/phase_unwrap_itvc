function [dx_map, dy_map] = calculate_gradients(map)
%CALCULATE_GRADIENTS Computes gradients of a 2D map (first order
%differences) x,y-directions correspond to column, row respectively.  
%   The derivative maps returned are of sizes Hx(W-1), (H-1)xW respectively
%   when compared to the original map's size (HxW)
%   dx_map(i,j) = map(i,j+1) - map(i,j)
%   dy_map(i,j) = map(i+1,j) - map(i,j)
dx_map = map(:,2:end) - map(:,1:end-1);
dy_map = map(2:end,:) - map(1:end-1,:);
end

