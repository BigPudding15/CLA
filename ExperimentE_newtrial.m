e = 10^(-4);

%% Dense case
n = 100;
A = hilb(n);
cond(A)

% Computing ground truth
x_gt = randn(n, 1);
b = A * x_gt;


% Solve system using Gauss_LS_pd
[x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, n_Gausspd] = Gauss_pd(A, b, e, 15, true);
% Solve system using CD_LS_pd
[x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, n_LSpd] = CD_LS_pd(A, b, e, 15, true);
mynGauss = n_Gausspd;
mynCD = n_LSpd;

mean_log_resid_G_LS_dense = log_resid_G_LS_dense;
mean_log_time_G_LS_dense  = log_time_G_LS_dense;
mean_log_flops_G_LS_dense = log_flops_G_LS_dense;

mean_log_resid_CD_LS_dense = log_resid_CD_LS_dense;
mean_log_time_CD_LS_dense = log_time_CD_LS_dense;
mean_log_flops_CD_LS_dense= log_flops_CD_LS_dense;

k = 5;

for j=1:k
    % Solve system using Gauss_LS_pd
    [x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, n_Gausspd] = Gauss_pd(A, b, e, 15, true);
    % Solve system using CD_LS_pd
    [x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, n_LSpd] = CD_LS_pd(A, b, e, 15, true);
    mynGauss = min(n_Gausspd,mynGauss);
    mynCD = min(n_LSpd,mynCD);
    mean_log_resid_G_LS_dense =    (1/(j+1))*(j*mean_log_resid_G_LS_dense(1:mynGauss)+log_resid_G_LS_dense(1:mynGauss));
    mean_log_time_G_LS_dense  = (1/(j+1))*(j*mean_log_time_G_LS_dense(1:mynGauss)+log_time_G_LS_dense(1:mynGauss));
    mean_log_flops_G_LS_dense = (1/(j+1))*(j*mean_log_flops_G_LS_dense(1:mynGauss)+log_flops_G_LS_dense(1:mynGauss));
    
    mynCD = min(mynCD,n_LSpd);
    mean_log_resid_CD_LS_dense = (1/(j+1))*(j*mean_log_resid_CD_LS_dense(1:mynCD)+log_resid_CD_LS_dense(1:mynCD));
    mean_log_time_CD_LS_dense = (1/(j+1))*(j*mean_log_time_CD_LS_dense(1:mynCD)+log_time_CD_LS_dense(1:mynCD));
    mean_log_flops_CD_LS_dense= (1/(j+1))*(j*mean_log_flops_CD_LS_dense(1:mynCD)+log_flops_CD_LS_dense(1:mynCD));
    
end

% Plot time - dense
subplot(1, 2, 1)
plot(mean_log_time_G_LS_dense, mean_log_resid_G_LS_dense);
hold on
plot(mean_log_time_CD_LS_dense, mean_log_resid_CD_LS_dense);
xlabel("time (s)")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
% Plot flops - dense
subplot(1, 2, 2)
plot(mean_log_flops_G_LS_dense, mean_log_resid_G_LS_dense);
hold on
plot(mean_log_flops_CD_LS_dense, mean_log_resid_CD_LS_dense);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off





%% Sparse case
n = 1000;
density = 1/log(n^2);
rc = 1/n;
A = sprandsym(n, density, rc, 1);
%% Computing ground truth
x_gt = randn(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS
[x_G_LS_sparse, log_resid_G_LS_sparse, log_x_G_LS_sparse, log_time_G_LS_sparse, log_flops_G_LS_sparse, ~] = Gauss_LS(A, b, e, 15, true);
%% Solve system using CD_LS
[x_CD_LS_sparse, log_resid_CD_LS_sparse, log_x_CD_LS_sparse, log_time_CD_LS_sparse, log_flops_CD_LS_sparse, ~] = CD_LS(A, b, e, 15, true);



%% Plot time - dense
subplot(2, 2, 1)
plot(log_time_G_LS_dense, log_resid_G_LS_dense);
hold on
plot(log_time_CD_LS_dense, log_resid_CD_LS_dense);
xlabel("time (s)")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - dense
subplot(2, 2, 2)
plot(log_flops_G_LS_dense, log_resid_G_LS_dense);
hold on
plot(log_flops_CD_LS_dense, log_resid_CD_LS_dense);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot time - sparse
subplot(2, 2, 3)
plot(log_time_G_LS_sparse, log_resid_G_LS_sparse);
hold on
plot(log_time_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("time")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - sparse
subplot(2, 2, 4)
plot(log_flops_G_LS_sparse, log_resid_G_LS_sparse);
hold on
plot(log_flops_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off