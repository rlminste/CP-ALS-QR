%% Figs3_4_sinsums
%
% Generates Figs 3 and 4 from
%   'CP Decomposition for Tensors via Alternating 
%    Least Squares with QR Decomposition'
%       - Minster, Viviano, Liu, Ballard
%
% tests accuracy of CP-ALS, CP-ALS-PINV, CP-ALS-QR, & CP-ALS-QR-SVD
% on 4 and 5-way sin of sums tensors

clear

%% 4-way (Fig 3)
d = 4; %number of modes
r = 4;
n = [64,128]; %dimension of tensor
maxiter = 40;
tol = 0;
rng(0)

figure,
for i = 1:2
    % form sin of sums tensor
    T = sinsum_full(d,n(i)); % rank 2^(d-1) representation

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
    subplot(1,2,i)
    semilogy(1:out_als.iters,out_als.relerr,'linewidth',2), hold on
    semilogy(1:out_pinv.iters,out_pinv.relerr,'--','linewidth',2)
    semilogy(1:out_qr.iters,out_qr.relerr,'-.','linewidth',2)
    semilogy(1:out_svd.iters,out_svd.relerr,':','linewidth',2)
    ylim([1e-12 10])
    legend('CP-ALS','CP-ALS-PINV','CP-ALS-QR','CP-ALS-QR-SVD')
    xlabel('iteration number')
    ylabel('relative error')
    title(sprintf('N = 4, n=%d',n(i))) 
    set(gca,'fontsize',16)
end

%% 5-way (Fig 4)
d = 5; %number of modes
r = 5;
n = [32,64]; %dimension of tensor
maxiter = 40;
tol = 0;
rng(0)

figure,
for i = 1:2
    % form sin of sums tensor
    T = sinsum_full(d,n(i)); % rank 2^(d-1) representation

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
    subplot(1,2,i)
    semilogy(1:out_als.iters,out_als.relerr,'linewidth',2), hold on
    semilogy(1:out_pinv.iters,out_pinv.relerr,'--','linewidth',2)
    semilogy(1:out_qr.iters,out_qr.relerr,'-.','linewidth',2)
    semilogy(1:out_svd.iters,out_svd.relerr,':','linewidth',2)
    ylim([1e-12 10])
    legend('CP-ALS','CP-ALS-PINV','CP-ALS-QR','CP-ALS-QR-SVD')
    xlabel('iteration number')
    ylabel('relative error')
    title(sprintf('N = 5, n=%d',n(i))) 
    set(gca,'fontsize',16)
end