function [fx,fy,fxx,fyx,fxy,fyy] = f_get_derivatives(f)
%F_GET_DERIVATIVES returns 1st and 2nd order derivatives of a 2D map
%   x-derivative along columns, y-derivative along rows. Simple forward
%   difference operator used
addpath ./lib/
[fx,fy] = calculate_gradients(f);
[fxx,fyx] = calculate_gradients(fx);
[fxy,fyy] = calculate_gradients(fy);
end

