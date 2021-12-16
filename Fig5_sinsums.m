%% Fig5_sinsums
%
% Generates Fig 5 from
%   'CP Decomposition for Tensors via Alternating 
%    Least Squares with QR Decomposition'
%       - Minster, Viviano, Liu, Ballard
%
% tests accuracy of CP-ALS, CP-ALS-PINV, CP-ALS-QR, & CP-ALS-QR-SVD
% on 10-way sin of sums tensor

clear

%% 10-way (left)
d = 10; %number of modes
r = 10;
n = 8; %dimension of tensor
maxiter = 40;
tol = 0;
rng(3)

% form sin of sums tensor
T = sinsums(d,n); % rank d representation

% Compute CP decomposition
% CP-ALS
[M_als,U_als,out_als1] = cp_als_time(T,r,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

% CP-ALS-PINV 
[M_pinv,U_pinv,out_pinv1] = cp_als_pinv(T,r,'init',U_als,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

% CP-ALS-QR
[M_qr,U_qr,out_qr1] = cp_als_qr(T,r,'init',U_als,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

% CP-ALS-QR-SVD
[M_svd,U_svd,out_svd1] = cp_als_qr_svd(T,r,'init',U_als,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');


%% 10-way (right)
d = 10; %number of modes
r = 10;
n = 8; %dimension of tensor
maxiter = 40;
tol = 0;
rng(0)

% form sin of sums tensor
T = sinsums(d,n); % rank d representation

% Compute CP decomposition
% CP-ALS
[M_als,U_als,out_als] = cp_als_time(T,r,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

% CP-ALS-PINV 
[M_pinv,U_pinv,out_pinv] = cp_als_pinv(T,r,'init',U_als,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

% CP-ALS-QR
[M_qr,U_qr,out_qr] = cp_als_qr(T,r,'init',U_als,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

% CP-ALS-QR-SVD
[M_svd,U_svd,out_svd] = cp_als_qr_svd(T,r,'init',U_als,'maxiters',maxiter,'tol',tol,'printitn',10,'errmethod','lowmem');

%% plot iterations vs relative error
figure,
subplot(1,2,1)
semilogy(1:out_als1.iters,out_als1.relerr,'linewidth',2), hold on
semilogy(1:out_pinv1.iters,out_pinv1.relerr,'--','linewidth',2)
semilogy(1:out_qr1.iters,out_qr1.relerr,'-.','linewidth',2)
semilogy(1:out_svd1.iters,out_svd1.relerr,':','linewidth',2)
ylim([1e-12 10])
legend('CP-ALS','CP-ALS-PINV','CP-ALS-QR','CP-ALS-QR-SVD')
xlabel('iteration number')
ylabel('relative error')
title(sprintf('N = 10, 1st initialization')) 
set(gca,'fontsize',16)

subplot(1,2,2)
semilogy(1:out_als.iters,out_als.relerr,'linewidth',2), hold on
semilogy(1:out_pinv.iters,out_pinv.relerr,'--','linewidth',2)
semilogy(1:out_qr.iters,out_qr.relerr,'-.','linewidth',2)
semilogy(1:out_svd.iters,out_svd.relerr,':','linewidth',2)
ylim([1e-12 10])
legend('CP-ALS','CP-ALS-PINV','CP-ALS-QR','CP-ALS-QR-SVD')
xlabel('iteration number')
ylabel('relative error')
title(sprintf('N = 10, 2nd initialization')) 
set(gca,'fontsize',16)