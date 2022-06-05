m = 1000;
n = 500; 
e = 10^(-4);


%% Dense case
A = randn(m, n);
%% Computing ground truth
x_gt = randn(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS
[x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, ~] = Gauss_LS(A, b, e, 300, true);
%% Solve system using GK
[x_GK_dense, log_resid_GK_dense, log_x_GK_dense, log_time_GK_dense, log_flops_GK_dense, ~] = GK(A, b, e, 300, true);
%% Solve system using RK
[x_RK_dense, log_resid_RK_dense, log_x_RK_dense, log_time_RK_dense, log_flops_RK_dense, ~] = RK(A, b, e, 300, true);
%% Solve system using CD_LS
[x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, ~] = CD_LS(A, b, e, 300, true);




%% Sparse case
density = 1 / log(m*n);
rc = 1 / (m*n)^(0.5);
A = sprandn(m, n, density,rc);
%% Computing ground truth
x_gt = rand(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS
[x_G_LS_sparse, log_resid_G_LS_sparse, log_x_G_LS_sparse, log_time_G_LS_sparse, log_flops_G_LS_sparse, ~] = Gauss_LS(A, b, e, 300, true);
%% Solve system using GK
[x_GK_sparse, log_resid_GK_sparse, log_x_GK_sparse, log_time_GK_sparse, log_flops_GK_sparse, ~] = GK(A, b, e, 300, true);
%% Solve system using RK
[x_RK_sparse, log_resid_RK_sparse, log_x_RK_sparse, log_time_RK_sparse, log_flops_RK_sparse, ~] = RK(A, b, e, 300, true);
%% Solve system using CD_LS
[x_CD_LS_sparse, log_resid_CD_LS_sparse, log_x_CD_LS_sparse, log_time_CD_LS_sparse, log_flops_CD_LS_sparse, ~] = CD_LS(A, b, e, 300, true);



%% Plot time - dense
subplot(2, 2, 1)
plot(log_time_G_LS_dense, log_resid_G_LS_dense);
hold on
plot(log_time_GK_dense, log_resid_GK_dense);
hold on
plot(log_time_RK_dense, log_resid_RK_dense);
hold on
plot(log_time_CD_LS_dense, log_resid_CD_LS_dense);
xlabel("time")
ylabel("error")
legend("Gauss LS", "GK", "RK", "CD LS")
hold off
%% Plot flops - dense
subplot(2, 2, 2)
plot(log_flops_G_LS_dense, log_resid_G_LS_dense);
hold on
plot(log_flops_GK_dense, log_resid_GK_dense);
hold on
plot(log_flops_RK_dense, log_resid_RK_dense);
hold on
plot(log_flops_CD_LS_dense, log_resid_CD_LS_dense);
xlabel("flops")
ylabel("error")
legend("Gauss LS", "GK", "RK", "CD LS")
hold off
%% Plot time - sparse
subplot(2, 2, 3)
plot(log_time_G_LS_sparse, log_resid_G_LS_sparse);
hold on
plot(log_time_GK_sparse, log_resid_GK_sparse);
hold on
plot(log_time_RK_sparse, log_resid_RK_sparse);
hold on
plot(log_time_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("time")
ylabel("error")
legend("Gauss LS", "GK", "RK", "CD LS")
hold off
%% Plot flops - sparse
subplot(2, 2, 4)
plot(log_flops_G_LS_sparse, log_resid_G_LS_sparse);
hold on
plot(log_flops_GK_sparse, log_resid_GK_sparse);
hold on
plot(log_flops_RK_sparse, log_resid_RK_sparse);
hold on
plot(log_flops_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("flops")
ylabel("error")
legend("Gauss LS", "GK", "RK", "CD LS")
hold off
