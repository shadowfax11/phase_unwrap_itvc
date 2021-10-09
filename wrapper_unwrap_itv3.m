function wrapper_unwrap_itv3(P,Pn,method,save_str)
%WRAPPER_UNWRAP_ITV3 Summary of this function goes here
%   Detailed explanation goes here
result_dir = './results/';
if strcmp(method,'IHTV')
    num_iter    = 5;
    eps1        = 1;
    eps2        = 0;
    lambda_1    = 1e0;
    lambda_2    = 0e0;
    beta_1      = 1e3;
    beta_2      = 0e3;
    mu_1        = 0e0;
    mu_2        = 0e0;
end
if strcmp(method,'ITVC')
    num_iter    = 5;
    eps1        = 1;
    eps2        = 1;
    lambda_1    = 1e0;
    lambda_2    = 1e0;
    beta_1      = 1e3;
    beta_2      = 1e3;
    mu_1        = 1e0;
    mu_2        = 1e0;
end
if strcmp(method,'ITVC-Base')
    num_iter    = 5;
    eps1        = 1;
    eps2        = 1;
    lambda_1    = 0e0;
    lambda_2    = 0e0;
    beta_1      = 0e3;
    beta_2      = 0e3;
    mu_1        = 1e0;
    mu_2        = 1e0;
end
if strcmp(method,'ITVC-TV')
    num_iter    = 5;
    eps1        = 1;
    eps2        = 1;
    lambda_1    = 1e0;
    lambda_2    = 1e0;
    beta_1      = 0e3;
    beta_2      = 0e3;
    mu_1        = 1e0;
    mu_2        = 1e0;
end
if strcmp(method,'ITVC-IRR')
    num_iter    = 5;
    eps1        = 1;
    eps2        = 1;
    lambda_1    = 0e0;
    lambda_2    = 0e0;
    beta_1      = 1e3;
    beta_2      = 1e3;
    mu_1        = 1e0;
    mu_2        = 1e0;
end
W_Pn = wrapToPi(Pn);
% preprocessing and preparing variables
sz = size(Pn);
P_est = zeros(sz);
res = wrapToPi(W_Pn - P_est);
est_i = {};
res_i = {};
init_val = P(1,1);
% quantities to track
res_mae = {};
res_rmse = {};
diff_mae = {};
diff_rmse = {};
num_2pi_jumps = {};

break_flag = 0;
for i=1:num_iter
    % perform one iteration of unwrapping
    [f_est, ~, ~] = unwrap_itv3(res,lambda_1,lambda_2,...
        beta_1,beta_2,mu_1,mu_2,eps1,eps2,init_val);
    x = 1:sz(1); y = 1:sz(2); [X,Y] = meshgrid(x,y);
    X(isnan(f_est)) = []; Y(isnan(f_est)) = []; f_est(isnan(f_est)) = [];
    x = 1:sz(1); y = 1:sz(2); [Xq,Yq] = meshgrid(x,y);
    f_est = reshape(griddata(X,Y,f_est,Xq(:),Yq(:)),sz(1),sz(2));
    P_est = P_est + f_est;
    res = wrapToPi(W_Pn - P_est);
    init_val = 0;
    % display results
    mae = mean(abs(res),'all');
    rmse = rms(res,'all');
    P_err = P - P_est;
    fprintf("%d out of %d iterations done\n",i,num_iter);
    fprintf("Residual stats \t MAE %.5f \t RMSE %.5f\n",mae,rmse);
    fprintf("Error stats \t MAE %.5f \t RMSE %.5f\n",mean(abs(P_err),'all'),rms(P_err,'all'));
    est_i{i} = P_est;
    res_i{i} = res;
    % track objective function costs
    res_mae{i} = mae;
    res_rmse{i} = rmse;
    diff_mae{i} = mean(abs(P_err),'all');
    diff_rmse{i} = rms(P_err,'all');
    num_2pi_jumps{i} = frac_2pi_jumps(res);
    if num_2pi_jumps{i} < 0.02
        fprintf("Stopping\n");
        break_flag = 1;
        break
    end
end
IMIN = 0;
if break_flag
    [fsim_val,~] = FeatureSIM(P,est_i{end});
    [mlv_val,~] = MLVSharpnessMeasure(est_i{end});
else
    FMIN = 1; 
    for j=1:numel(num_2pi_jumps)
        if FMIN > num_2pi_jumps{j}
            FMIN = num_2pi_jumps{j};
            IMIN = j;
        end
    end
    [fsim_val,~] = FeatureSIM(P,est_i{IMIN});
    [mlv_val,~] = MLVSharpnessMeasure(est_i{IMIN});
end
fprintf("Metrics: \n MAE  %.5f \t RMSE %.5f \n FSIM %.5f \n MLV  %.5f\n",...
    mae,rmse,fsim_val,mlv_val);

save([result_dir,save_str,'-',method,'.mat'],'est_i','res_i','P','Pn',...
    'res_mae','res_rmse','diff_mae','diff_rmse','IMIN',...
    'lambda_1','lambda_2','beta_1','beta_2','mu_1','mu_2','fsim_val','mlv_val');

end

function f = frac_2pi_jumps(res)
[dx,dy] = calculate_gradients(res);
Nx = dx~=wrapToPi(dx);
Ny = dy~=wrapToPi(dy);
Nxy = zeros(size(res));
Nxy(1:end-1,1:end-1) = Nx(1:end-1,:)|Ny(:,1:end-1);
Nxy(end,1:end-1) = Nx(end,:);
Nxy(1:end-1,end) = Ny(:,end);
f = sum(Nxy,'all')/(numel(res));
end
