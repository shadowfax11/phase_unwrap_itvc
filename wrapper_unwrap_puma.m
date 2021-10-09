function wrapper_unwrap_puma(P,Pn,eps,method_name,save_str)
%WRAPPER_UNWRAP_PUMA Summary of this function goes here
%   Detailed explanation goes here
addpath ./lib/
method = method_name;
result_dir = '../results_puma/';
num_iter = 1;
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
    [f_est,~,~] = puma_ho(res,eps);
    if i==1
        f_est = f_est - f_est(1,1) + P(1,1); 
    else
        f_est = f_est - f_est(1,1) + res(1,1);
    end
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
    if num_2pi_jumps{i} < 0.000
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
    'eps','fsim_val','mlv_val');
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
