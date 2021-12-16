%% Tab2_collinear
%
% Generates figures for Table 2 from
%   'CP Decomposition for Tensors via Alternating 
%    Least Squares with QR Decomposition'
%       - Minster, Viviano, Liu, Ballard
%
% tests accuracy of CP-ALS, CP-ALS-PINV, CP-ALS-QR, & CP-ALS-QR-SVD
% on random tensors with collinear factor matrices and added noise

% parameters
n = 50; % dimension
sz = [n;n;n]; % 3-mode tensor
maxiter = 200;
tol = 0;
ntrials = 10;
noise = [1e-4,1e-7,1e-10]; % noise
collinear = [1-1e-4,1-1e-7,1-1e-10]; % collinearity of factor matrices
r = 5; % rank

scores_all = cell(3,3);
iters_all = cell(3,3);
relerr_all = cell(3,3);

for i = 1:3
    ns = noise(i);
    for j = 1:3
        coll = collinear(j);
        
        % create ktensor with collinear factor matrices
        info = create_problem('Size',sz,'Num_Factors',r,'Noise',ns,'Factor_Generator',@(m,n) matrandcong(m,n,coll));
        T = info.Data; % tensor
        Mtrue = info.Soln; % exact answer (ktensor)
        
        scores = zeros(ntrials,5);
        iters = zeros(ntrials,5);
        relerr = zeros(ntrials,5);
        
        for k = 1:ntrials
            % CP with Gauss-Newton
            [err_gn,scores_gn,iters_gn] = cp_gn(T,r,Mtrue);
            scores(k,1) = scores_gn;
            iters(k,1) = iters_gn;
            relerr(k,1) = err_gn;

            % CP-ALS
            [M_als,U_als,out_als] = cp_als_times(T,r,'maxiters',maxiter,'printitn',100,'tol',tol,'errmethod','full');
            scores(k,2) = score(M_als,Mtrue,'lambda_penalty',false,'greedy',false);
            iters(k,2) = out_als.iters;
            relerr(k,2) = out_als.relerr(out_als.iters);

            % CP-ALS-PINV
            [M_pinv,~,out_pinv] = cp_als_pinv(T,r,'init',U_als,'maxiters',maxiter,'printitn',100,'tol',tol,'errmethod','full');
            scores(k,3) = score(M_pinv,Mtrue,'lambda_penalty',false,'greedy',false);
            iters(k,3) = out_pinv.iters;
            relerr(k,3) = out_pinv.relerr(out_pinv.iters);

            % CP-ALS-QR
            [M_qr,~,out_qr] = cp_als_qr(T,r,'init',U_als,'maxiters',maxiter,'printitn',100,'tol',tol,'errmethod','full');
            scores(k,4) = score(M_qr,Mtrue,'lambda_penalty',false,'greedy',false);
            iters(k,4) = out_qr.iters;
            relerr(k,4) = out_qr.relerr(out_qr.iters);

            % CP-ALS-QR-SVD
            [M_svd,~,out_svd] = cp_als_qr_svd(T,r,'init',U_als,'maxiters',maxiter,'printitn',100,'tol',tol,'errmethod','full');
            scores(k,5) = score(M_svd,Mtrue,'lambda_penalty',false,'greedy',false);
            iters(k,5) = out_svd.iters;
            relerr(k,5) = out_svd.relerr(out_svd.iters);
        end
        
        scores_all{i,j} = scores;
        iters_all{i,j} = iters;
        relerr_all{i,j} = relerr;
    end
end
      
%% make boxplots
for i = 1:3
    for j = 1:3
        figure,
        subplot(3,1,1)
        boxplot(scores_all{i,j},'Labels',{'','','','',''})
        ylabel('Scores')
        ylim([-.1 1])
        title(sprintf('noise = %d, coll = %d',noise(i),collinear(j)))
        set(gca,'fontsize',16)
        
        subplot(3,1,2)
        boxplot(iters_all{i,j},'Labels',{'','','','',''})
        ylabel('Iterations')
        ylim([-10 510])
        set(gca,'fontsize',16)

        subplot(3,1,3)
        boxplot(relerr_all{i,j},'Labels',{'Gauss-Newton','CP-ALS','CP-ALS-PINV','CP-ALS-QR','CP-ALS-QR-SVD'})
        ylabel('Relative Error')
        ylim([1e-11 1e-2])
        set(gca,'fontsize',16)
    end
end