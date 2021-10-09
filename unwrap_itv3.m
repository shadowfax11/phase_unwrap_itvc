

function [map,residual,costs] = unwrap_itv3(wrapped_map,lambda_1,lambda_2,...
    beta_1,beta_2,mu_1,mu_2,eps1,eps2,init_val)
%UNWRAP_ITV3 Performs iTV phase unwrapping (combined optimization of
%level-1 and level-2) for one iteration
%   Uses admm_solver_v5() function 
[dx, dy] = calculate_gradients(wrapped_map);
[dxx,dyx] = calculate_gradients(dx);
[dxy,dyy] = calculate_gradients(dy);
g1 = wrapToPi(dx);
g2 = wrapToPi(dy);
h1 = wrapToPi(dxx); 
h2 = wrapToPi(dyx);
h3 = wrapToPi(dxy);
h4 = wrapToPi(dyy);
[fx,fy,costs] = admm_solver_v5(g1, g2, h1, h2, h3, h4, ...
    [lambda_1,lambda_2,beta_1,beta_2,mu_1,mu_2,eps1,eps2]);
map = integrate_derivatives(fx,fy,2,init_val);
% map = integrate_derivatives(fx,fy,1,init_val);
residual = wrapToPi(wrapped_map - map);
end

