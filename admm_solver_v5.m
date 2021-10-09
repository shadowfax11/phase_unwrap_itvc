function [u1,u2,costs] = admm_solver_v5(g1,g2,h1,h2,h3,h4,hyperparam_lst)
%ADMM_SOLVER_V5 This function sets up the ADMM solver for minimizing -
%min_u,v ||u-g||_2^2 + ||v-h||_2^2 + lambda_1*HTV(u) + lambda_2*HTV(v) +
%lambda1*beta_1*||Mu||_1^1 | + lambda2*beta_2*||Mv_12||_1^1 + lambda2*beta_2*||Mv_34||_1^1 +
%mu_1*||Ev_23||_2^2 + mu_2*||Du - v||_2^2
%   HTV operator is a 2nd order derivative operator (L-1 norm)
%   M is the irrotationality operator
%   E is the equality operator
%   D is the forward difference operator
%   g = [g1;g2]     h = [h1;h2;h3;h4]
%   u = [u1;u2]     v = [v1;v2;v3;v4]
%   v_12 = [v1;v2]  v_23 = [v2;v3]  v_34 = [v3;v4]

% preparing hyperparameters
lambda_1 = hyperparam_lst(1);
lambda_2 = hyperparam_lst(2);
beta_1 = hyperparam_lst(3);
beta_2 = hyperparam_lst(4);
mu_1 = hyperparam_lst(5);
mu_2 = hyperparam_lst(6);
eps1 = hyperparam_lst(7); 
eps2 = hyperparam_lst(8);

% preprocessing and preparing variables
s_g1 = size(g1); s_g2 = size(g2); 
s_h1 = size(h1); s_h2 = size(h2); s_h3 = size(h3); s_h4 = size(h4);

n_g1 = numel(g1); n_g2 = numel(g2); 
n_h1 = numel(h1); n_h2 = numel(h2); n_h3 = numel(h3); n_h4 = numel(h4);

g1 = reshape(g1, [n_g1,1]);
g2 = reshape(g2, [n_g2,1]);
h1 = reshape(h1, [n_h1,1]);
h2 = reshape(h2, [n_h2,1]);
h3 = reshape(h3, [n_h3,1]);
h4 = reshape(h4, [n_h4,1]);

% create HTV operators
F_g12 = create_htv_operators(s_g1, s_g2, n_g1, n_g2);
F_h12 = create_grad_operator(s_h1, s_h2, n_h1, n_h2);
F_h34 = create_grad_operator(s_h3, s_h4, n_h3, n_h4);
n_htv_g12 = size(F_g12,1);
n_htv_h12 = size(F_h12,1);
n_htv_h34 = size(F_h34,1);

% create IRR operators
F_g12 = [F_g12; (beta_1/lambda_1)*create_irr_operator(s_g1, s_g2, n_g1, n_g2)]; 
F_h12 = [F_h12; (beta_2/lambda_2)*create_irr_operator(s_h1, s_h2, n_h1, n_h2)];
F_h34 = [F_h34; (beta_2/lambda_2)*create_irr_operator(s_h3, s_h4, n_h3, n_h4)];

S_g12 = size(F_g12); 
S_h12 = size(F_h12);
S_h34 = size(F_h34);

F_g12 = [F_g12, sparse(S_g12(1),S_h12(2)+S_h34(2))];
F_h12 = [sparse(S_h12(1),S_g12(2)), F_h12, sparse(S_h12(1),S_h34(2))];
F_h34 = [sparse(S_h34(1),S_g12(2)+S_h12(2)), F_h34];

% create A matrix components
A1 = sparse(1:n_g1+n_g2+n_h1+n_h2+n_h3+n_h4,1:n_g1+n_g2+n_h1+n_h2+n_h3+n_h4,...
    [eps1*ones(1,n_g1+n_g2),eps2*ones(1,n_h1+n_h2+n_h3+n_h4)],...
    n_g1+n_g2+n_h1+n_h2+n_h3+n_h4,n_g1+n_g2+n_h1+n_h2+n_h3+n_h4);
i = [1:n_h2, 1:n_h2];
j = n_g1+n_g2+n_h1+[1:n_h2, n_h2+1:n_h2+n_h3];
v = sqrt(mu_1)*[ones(1,n_h2), -1*ones(1,n_h2)];
A2e = sparse(i,j,v,n_h2,n_g1+n_g2+n_h1+n_h2+n_h3+n_h4,2*n_h2);
D_u = create_grad_operator(s_g1,s_g2,n_g1,n_g2);
D_u = sqrt(mu_2)*[D_u, -1*speye(n_h1+n_h2+n_h3+n_h4)];
A2 = [A2e; D_u];

% create b-vector
b = [g1;g2;h1;h2;h3;h4;zeros(size(A2,1),1)];

% perform ADMM 
[u_hat, ~] = f_solver_admm_v5(A1, A2, b, lambda_1, lambda_2, F_g12, F_h12, F_h34, 1);

% print costs of each term
% fprintf("Lvl1 L-2 cost: %.5f\n",...
costs(1) =    eps1*0.5*norm(u_hat(1:n_g1+n_g2)-b(1:n_g1+n_g2))^2;
%     );
% fprintf("Lvl2 L-2 cost: %.5f\n",...
costs(2) =    eps2*0.5*norm(u_hat(n_g1+n_g2+1:n_g1+n_g2+n_h1+n_h2+n_h3+n_h4)-b(n_g1+n_g2+1:n_g1+n_g2+n_h1+n_h2+n_h3+n_h4))^2;
% fprintf("EqlM L-2 cost: %.5f\n",...
costs(3) =    0.5*norm(A2e*u_hat)^2;
% fprintf("Grad L-2 cost: %.5f\n",...
costs(4) =    0.5*norm(D_u*u_hat)^2;
% fprintf("Lvl1 HTV cost: %.5f\n",...
costs(5) =    lambda_1*norm(F_g12(1:n_htv_g12,:)*u_hat,1);
% fprintf("Lvl2 HTV cost: %.5f\n",...
%costs(6) =    lambda_2*(norm((F_h12(1:n_htv_h12,:)*u_hat) + (F_h34(1:n_htv_h34,:)*u_hat),1));
% fprintf("Lvl1 Irr cost: %.5f\n",...
%costs(7) =    lambda_1*norm(F_g12(n_htv_g12+1:end,:)*u_hat,1);
% fprintf("Lvl2 Irr cost: %.5f\n",...
%costs(8) =    lambda_2*(norm((F_h12(n_htv_h12+1:end,:)*u_hat) + (F_h34(n_htv_h34+1:end,:)*u_hat),1));

u1 = u_hat(1:n_g1);
u2 = u_hat(n_g1+1:n_g1+n_g2);
u1 = reshape(u1, s_g1);
u2 = reshape(u2, s_g2);
end

function D = create_grad_operator(s1,s2,n1,n2)
k1 = s1(1);
i = [1:n1, n1+1:2*n1, 1:n1-1, n1+1:2*n1-k1];
j = [1:n1, 1:n1, 2:n1, 1+k1:n1];
v = [-1*ones(1,2*n1), ones(1,n1-1), ones(1,n1-k1)];
D1 = sparse(i,j,v,2*n1,n1);
D1(k1:k1:n1,:) = 0;
D1(2*n1-k1+1:end,:) = 0;
sz = size(D1);
D1 = [D1((sz(1)/2)+1:end,:); D1(1:sz(1)/2,:)];
k2 = s2(1);
i = [1:n2, n2+1:2*n2, 1:n2-1, n2+1:2*n2-k2];
j = [1:n2, 1:n2, 2:n2, 1+k2:n2];
v = [-1*ones(1,2*n2), ones(1,n2-1), ones(1,n2-k2)];
D2 = sparse(i,j,v,2*n2,n2);
D2(k2:k2:n2,:) = 0;
D2(2*n2-k2+1:end,:) = 0;
sz = size(D2);
D2 = [D2((sz(1)/2)+1:end,:); D2(1:sz(1)/2,:)];
z0 = sparse(size(D1,1),n2);
z3 = sparse(size(D2,1),n1);
D = [D1, z0; z3, D2];
D = D(any(D,2),:);
end

function D = create_htv_operators(s1,s2,n1,n2)
k1 = s1(1);
i = [1:n1, n1+1:2*n1, 1:n1-1, n1+1:2*n1-k1];
j = [1:n1, 1:n1, 2:n1, 1+k1:n1];
v = [-1*ones(1,2*n1), ones(1,n1-1), ones(1,n1-k1)];
D1 = sparse(i,j,v,2*n1,n1);
D1(k1:k1:n1,:) = 0;
D1(2*n1-k1+1:end,:) = 0;
z0 = sparse(2*n1,n1);
z1 = sparse(4*n1,n2);
z2 = sparse(4*n2,n1);
z3 = sparse(2*n2,n2);
D1 = [D1, z0;z0, D1]*D1;
k2 = s2(1);
i = [1:n2, n2+1:2*n2, 1:n2-1, n2+1:2*n2-k2];
j = [1:n2, 1:n2, 2:n2, 1+k2:n2];
v = [-1*ones(1,2*n2), ones(1,n2-1), ones(1,n2-k2)];
D2 = sparse(i,j,v,2*n2,n2);
D2(k2:k2:n2,:) = 0;
D2(2*n2-k2+1:end,:) = 0;
D2 = [D2, z3; z3, D2]*D2;
D = [D1, z1; z2, D2];
end

function M = create_irr_operator(s1,s2,n1,n2)
k = s1(2)*s2(1);
i = [1:k,1:k,1:k,1:k];
j11 = 1:s2(1);
j1 = 1:s2(1);
for p=1:s1(2)-1
    j1 = [j1, p*s1(1)+j11];
end
j2 = s1(1)*s1(2) + s2(1) + [1:k];
j3 = j1+1;
j4 = s1(1)*s1(2) + [1:k];
j = [j1, j2, j3, j4];
v = [ones(1,2*k),-1*ones(1,2*k)];
M = sparse(i,j,v,k,n1+n2,4*k);
end

