function [x_hat,z_hat] = f_solver_admm_v5(A1,A2,b,lambda1,lambda2,F1,F2,F3,rho,~)
%F_SOLVER_ADMM This function uses ADMM to solve a minimization problem
%   

ABS_TOL = 1e-4;
REL_TOL = 1e-2;
MAX_ITERATIONS = 20;

m = size(A1,1);
n1 = size(F1,1);
n2 = size(F2,1);
n3 = size(F3,1);

% initialize solutions/variables
x = b(1:m);
b0 = b(1:m);
u1 = zeros(n1,1); u2 = zeros(n2,1); u3 = zeros(n3,1);
z1 = zeros(n1,1); z2 = zeros(n2,1); z3 = zeros(n3,1);

%preprocessing to ease calculations
% F = [F1;F2];
f_A1 = [];
f_A2 = [];
f_F1 = [];
f_F2 = [];
f_F3 = [];
% fh1 = figure; title("Costs");

for k=1:MAX_ITERATIONS
    
    % x-minimization update
    t = b0 + rho*((double(z1-u1))'*F1)' + rho*((double(z2-u2))'*F2)' + rho*((double(z3-u3))'*F3)';
    [x_new,~] = bicgstab(@afun,t,1e-3,200,[],[],x);
    
    % z-minimization update
    htv1 = F1*x_new;
    htv2 = F2*x_new;
    htv3 = F3*x_new;
    z_new1 = soft_threshold(htv1 + (u1/rho), lambda1/rho);
    z_new2 = soft_threshold(htv2 + (u2/rho), lambda2/rho);
    z_new3 = soft_threshold(htv3 + (u3/rho), lambda2/rho);
    % u-minimization update / y-minimization update
    res_p1 = htv1 - z_new1;
    res_p2 = htv2 - z_new2;
    res_p3 = htv3 - z_new3;
    u1 = u1 + rho*(res_p1);
    u2 = u2 + rho*(res_p2);
    u3 = u3 + rho*(res_p3);
    % error level checking
    primal_residue = norm([res_p1;res_p2;res_p3]);
    dual_residue = norm(-rho * (((z_new1-z1)'*F1)') -rho * (((z_new2-z2)'*F2)') -rho * (((z_new3-z3)'*F3)'));
    eps_primal = sqrt(n1+n2+n3)*ABS_TOL + REL_TOL*max(norm([htv1;htv2;htv3]), norm([-z_new1;-z_new2;-z_new3]));
    eps_dual = sqrt(m)*ABS_TOL + REL_TOL*norm(rho*((u1' * F1)')+rho*((u2' * F2)')+rho*((u3' * F3)'));
    
    x = x_new;
    z1 = z_new1;
    z2 = z_new2;
    z3 = z_new3;
    
    f_A1(k) = 0.5*norm(A1*x - b0, 'fro')^2;
    f_A2(k) = 0.5*norm(A2*x, 'fro')^2;
    f_F1(k) = lambda1*sum(abs(F1*x));
    f_F2(k) = lambda2*sum(abs(F2*x));
    f_F3(k) = lambda2*sum(abs(F3*x));
%     fprintf("%.5f %.5f %.5f\n",f_A1(k),f_A2(k),f_F1(k));
%     figure(fh1); 
%     semilogy(1:k,[f_A1', f_A2', f_F1', f_F2', f_F3']);
%     legend('A1','A2','F1','F2','F3'); title("Costs");
%     drawnow;
    if (primal_residue < eps_primal && ...
       dual_residue < eps_dual)
         break;
    end
%     fprintf("%d out of %d iterations done\n",k,MAX_ITERATIONS);
end

% final estimations
x_hat = x;
z_hat = [z1;z2;z3];

    function y = afun(x)
        % change it to make it faster
%         y = x + (rho*(F'*(F*x)));
        vI = A1*x;
        vA = A2*x;
        v1 = F1*x;
        v2 = F2*x;
        v3 = F3*x;
        y = (vI'*A1)' + (vA'*A2)' + rho*(v1'*F1)' + rho*(v2'*F2)' + rho*(v3'*F3)';
    end

end

% This function performs the soft-thresholding operation
function v = soft_threshold(v,thresh,w)
if nargin==3
    thresh = thresh*w;
end
v = max(0,v-thresh) - max(0,-v-thresh);
end
