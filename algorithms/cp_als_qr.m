function [P,Uinit,output] = cp_als_qr(X,R,varargin)
%CP_ALS_QR Compute a CP decomposition of any type of tensor.
%
%   M = CP_ALS(X,R) computes an estimate of the best rank-R
%   CP model of a tensor X using an alternating least-squares
%   algorithm.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor. The result P is a ktensor.
%
%   M = CP_ALS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations; 0 for no printing {1}
%
%   [M,U0] = CP_ALS(...) also returns the initial guess.
%
%   [M,U0,out] = CP_ALS(...) also returns additional output that contains
%   the input parameters.
%
%   Note: The "fit" is defined as 1 - norm(X-full(P))/norm(X) and is
%   loosely the proportion of the data described by the CP model, i.e., a
%   fit of 1 is perfect.
%
%   NOTE: Updated in various minor ways per work of Phan Anh Huy. See Anh
%   Huy Phan, Petr Tichavsk?, Andrzej Cichocki, On Fast Computation of
%   Gradients for CANDECOMP/PARAFAC Algorithms, arXiv:1204.1586, 2012.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   M = cp_als(X,2);
%   M = cp_als(X,2,'dimorder',[3 2 1]);
%   M = cp_als(X,2,'dimorder',[3 2 1],'init','nvecs');
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   [M,U0,out] = cp_als(X,2,'dimorder',[3 2 1],'init',U0);
%   M = cp_als(X,2,out.params); %<-- Same params as previous run
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.
%
% Adapted to use a QR to solve the LS problems instead of normal
% equations, also includes timing of iteration parts, and outputs relative
% error and fit. Used in 'CP Decomposition for Tensors via Alternating Least Squares with QR
% Decomposition'
% - Adapted by Irina Viviano & Rachel Minster, 2021




%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('tol',1e-4,@isscalar);
params.addParameter('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParameter('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.addParameter('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParameter('printitn',1,@isscalar);
params.addParameter('errmethod','fast',@(x) ismember(x,{'fast','full','lowmem'}));
params.parse(varargin{:});

%% Copy from params object
fitchangetol = params.Results.tol;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;
errmethod = params.Results.errmethod;


%% Error checking 

%% Set up and error checking on initial guess for U.
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end)
        if ~isequal(size(Uinit{n}),[size(X,n) R])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = rand(size(X,n),R);
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = nvecs(X,n,R);
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

if printitn>0
  fprintf('\nCP_ALS_QR:\n');
end


%% Main Loop: Iterate until convergence

%%% Changes for cp_als_qr start here: %%%

%%% Initialize a cell array Qs and Rs to hold decompositions of factor matrices. %%%
Qs = cell(N,1); %%% The Kronecker product of these tells us part of the Q of the Khatri-Rao product. %%%
Rs = cell(N,1); %%% The Khatri-Rao product of these tells us the rest of Q and the R of the Khatri-Rao product. %%%

%%% Compute economy-sized QR decomposition. %%%
for i = 1:N
    if ~isempty(U{i})
        [Qs{i}, Rs{i}] = qr(U{i},0); 
    end
end
   
for iter = 1:maxiters
    t_ttm = 0; % TTM
    t_qrf = 0; % QR of factor matrices
    t_kr = 0; % Computing Q0
    t_q0 = 0; % Applying Q0
    t_back = 0;
    t_lamb = 0;
    t_err = 0;
        
    fitold = fit;
        
    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)
       
        %%% Compute the QR of the Khatri-Rao product of Rs. %%%
        %%% First compute the Khatri-Rao product on all modes but n. %%%
        tic;  M = khatrirao(Rs{[1:n-1,n+1:N]},'r'); t = toc; t_kr = t_kr + t;
        
        %%% Compute the explicit QR factorization.
        tic;  [Q0,R0] = qr(M,0); t = toc; t_kr = t_kr + t;

        %%% TTM on all modes but mode n. %%%
        tic; Y = ttm(X,Qs,-n,'t'); t = toc; t_ttm = t_ttm + t;
        
        %%% Now multiply by Q0 on the right. %%%
        
        if isa(Y,'ktensor')
            %%% For a ktensor: %%%
            %%% Save all the factor matrices of Y in a cell array. %%%
            %%% Then, we can compute the Khatri Rao product in one line.%%%
            K = cell(N,1);
            
            for k = 1:N
                %if k ~= n
                    K{k} = Y.U{k};
                %end
            end
            
            %%% Apply Q0
            tic; Z = Y.U{n} * (khatrirao(K{[1:n-1,n+1:N]},'r')' * Q0); t = toc; t_q0 = t_q0 + t;

            %%% Calculate updated factor matrix by backsolving with R0' and Z. %%%
            tic; U{n} = double(Z) / R0'; t = toc; t_back = t_back + t;
            
        else
            %%% For any other tensor: %%%
            %%% Apply Q0 %%%
            tic; Z = tenmat(Y,n) * Q0; t = toc; t_q0 = t_q0 + t;

            %%% Calculate updated factor matrix by backsolving with R0' and Z. %%%
            tic; U{n} = double(Z) / R0'; t = toc; t_back = t_back + t;

        end
               
        % Normalize each vector to prevent singularities in coefmatrix
        tic;
        if iter == 1
            lambda = sqrt(sum(U{n}.^2,1))'; %2-norm
        else
            lambda = max( max(abs(U{n}),[],1), 1 )'; %max-norm
        end 
        
        Unew = bsxfun(@rdivide, U{n}, lambda'); t = toc; t_lamb = t_lamb + t;
        U{n} = Unew;
        
        %%% Recompute QR factorization for updated factor matrix. %%%
        tic; [Qs{n}, Rs{n}] = qr(U{n},0); t_qrf = toc;
    end

    %%% Changes for cp_als_qr end here. %%%
        
    P = ktensor(lambda,U);
    
    tic;
    if normX == 0
       Rscaled = lambda.*R0';
       prod = U{dimorder(end)}*Rscaled;
       iprod = Z(:)'*prod(:);
       fit = norm(P)^2 - 2 * iprod;
    else
       switch errmethod
            case 'fast'
                % fast inner product calculation
                Rscaled = R0.*lambda';
                prod = U{dimorder(end)}*Rscaled';
                iprod = Z(:)'*prod(:); 

                % fast norm(P) calculation: < Lambda R0^T R0 Lambda, Rs{N}^T Rs{N} >
                RscaledGram = Rscaled'*Rscaled;
                RnGram = Rs{dimorder(end)}'*Rs{dimorder(end)};
                normPsq = RscaledGram(:)'*RnGram(:);
                
                normresidual = sqrt( abs(normX^2 + normPsq - 2 * iprod) );
            case 'full'
                normresidual = norm(full(X) - full(P));
            case 'lowmem'
                normresidual = normdiff(X,P); 
        end
        fit = 1 - (normresidual / normX); %fraction explained by model
        %%% Change this to just be relative error to see the error go down. %%%
        rel_err(iter,:) = normresidual / normX; 
    end
    fitchange = abs(fitold - fit); t = toc; t_err = t_err + t;
        
    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        flag = 0;
    else
        flag = 1;
    end
        
    %%% If the fit is NaN, just stop the process. %%%
    if isnan(fit)
        break;
    end
    
    if (mod(iter,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
    end

    % Check for convergence
    if (flag == 0)
        break;
    end     
    
    times(iter,:) = [t_ttm, t_qrf, t_kr, t_q0, t_back, t_lamb, t_err];
end   


%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
P = fixsigns(P);



if printitn>0
    if normX == 0
        fit = norm(P)^2 - 2 * innerprod(X,P);
    else
        switch errmethod
            case 'fast'
                Rscaled = lambda.*R0';
                prod = U{dimorder(end)}*Rscaled;
                iprod = Z(:)'*prod(:);

                % fast norm(P) calculation: < Lambda R0^T R0 Lambda, Rs{N}^T Rs{N} >
                RscaledGram = Rscaled*Rscaled';
                RnGram = Rs{dimorder(end)}'*Rs{dimorder(end)};
                normPsq = RscaledGram(:)'*RnGram(:);
                normresidual = sqrt( abs(normX^2 + normPsq - 2 * iprod) );
            case 'full'
                normresidual = norm(full(X) - full(P));
            case 'lowmem'
                normresidual = normdiff(X,P);
        end
        fit = 1 - (normresidual / normX); %fraction explained by model
        rel_err(iter,:) = normresidual / normX;
    end
  fprintf(' Final f = %e \n', fit);
end



output = struct;
output.params = params.Results;
output.iters = iter;
output.relerr = rel_err; %%% Add a rel_err vector to output
output.fit = fit;
output.times = times;

