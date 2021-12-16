function [err,scores,iter] = cp_gn(T,r,Mtrue)
% wrapper to use GN on CP
% pass in true solution Mtrue to calculate score

addpath('tensor_toolbox-master')
addpath('tensorlab')

%% Configuring options
options.Display = true; % Show progress on the command line.
options.Compression = false;
options.Algorithm = @cpd_nls;
options.PlaneSearch = false;
options.Initialization = @cpd_rnd;
options.LineSearch = false;
options.Refinement = false;
options.MaxIter = 500;
options.AlgorithmOptions.TolFun = 1e-15; % Set function tolerance stop criterion
options.AlgorithmOptions.TolX = 1e-15; % Set step size tolerance stop criterion
% options.AlgorithmOptions.TolLargeScale = inf;
%% CP with Gauss Newton
[U_gn,output] = cpd(double(T), r, options);

% compute error, score, save number of iterations
err = output.Algorithm.relerr;
iter = output.Algorithm.iterations;
scores = score(arrange(ktensor(U_gn)),Mtrue,'lambda_penalty',false,'greedy',false);

end