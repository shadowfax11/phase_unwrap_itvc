function [curl_map] = calculate_curl(dx_map,dy_map)
%CALCULATE_CURL Given the gradient field maps, calculates the curl (2x2
%loop summation)
%   curl_map(i,j) = dx_map(i,j)+dy_map(i,j+1)-dx_map(i+1,j)-dy_map(i,j);
%   dx_map is of size Hx(W-1)
%   dy_map is of size (H-1)xW
%   curl_map is of size (H-1)x(W-1)
assert(size(dx_map,1)-1==size(dy_map,1));
assert(size(dx_map,2)==size(dy_map,2)-1);
curl_map = dx_map(1:end-1,:) + ...
            dy_map(:,2:end) - ...
            dx_map(2:end,:) - ...
            dy_map(:,1:end-1);
end

