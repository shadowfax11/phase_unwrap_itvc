function [map, dx, dy] = integrate_lvl2(dxx,dyx,dxy,dyy,init_val1,init_val2,init_val)
%INTEGRATE_LVL2 Summary of this function goes here
%   Detailed explanation goes here
dx = integrate_derivatives(dxx,dyx,2,init_val1);
dy = integrate_derivatives(dxy,dyy,2,init_val2);
map = integrate_derivatives(dx,dy,2,init_val);
end

