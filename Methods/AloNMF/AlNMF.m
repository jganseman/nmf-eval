function [W, H, HIS] = AlNMF(V, opts)
% Nonnegative Matrix Factorization (NMF) via Anti-lopsided + Greedy Coordinate Descent
% This software solve the following least squares NMF problem:
%
%	min_{W,H} 1/2*||V-WH||_F^2 + mu_1*||W||_1 + mu_2*||W||_2^2 + beta_1*||H||_1 + beta_2*||H||_2^2    s.t. W>=0, H>=0
%
% input: 
%		V: the input n by m dense matrix
%		k: the specified rank
%		opts.maxIter: maximum number of iterations
%		opts.W0: initial of W (n by k dense matrix)
%		opts.H0: initial of H (k by m dense matrix)
%		opts.verbose: 1: compute objective value per iteration.
%			   0: do not compute objective value per iteration. (default)
%
% output: 
%		NMF_GCD will output nonnegative matrices W, H, such that WH is an approximation of V
%		W: n by k dense matrix
%		H: k by m dense matrix
%		HIS.errors: objective values verus iterations. 
%		HIS.time: time verus iterations. 
%
maxIter = opts.maxIter;
Winit = opts.W0;
Hinit = opts.H0;

W = Winit;
H = Hinit;
k = size(Winit, 2);

%% Stopping tolerance for subproblems
tolW = 1e-2;
tolH = 1e-2;
maxIterW = k;
maxIterH = k;

startTime = tic;
HIS.time = toc(startTime);
HIS.error = norm(V - W*H, 'fro')^2/2.0;
[gradW, gradH] = getGrad(V, W, H);
normGradW = normGrad(gradW, W);
normGradH = normGrad(gradH, H);
HIS.gradW = normGradW;
HIS.gradH = normGradH;
HIS.iterW = 0.0;
HIS.iterH = 0.0;
HIS.iter = 0;
HIS.grad = (normGradW+normGradH);
fprintf('Begining: %.4E\n', norm(V-W*H, 'fro')^2/2.0);
W=W';

maxNumberThreads = opts.maxNumberThreads + 0.0;
[n, m] = size(V);

for i = 1:maxIter,
    startTime = tic;
    % update variables of H
    WW = W*W'; % (nxk)
    VW = -W*V; % (kxn) * (nxm)
    if isfield(opts, 'mu_1'), VW = VW + opts.mu_1; end
    if isfield(opts, 'mu_2'), WW = WW + 2 * opts.mu_2 * eye(k); end
    [H, iterH] = doAlIter(WW, VW, H, maxIterH, tolH, 0, maxNumberThreads);
    
    HH = H*H'; %(kxm)*(mxk)
    VH = -H*V'; %(kxm)*(mxn)]
    if isfield(opts, 'beta_1'),   VH = VH + opts.beta_1;  end
    if isfield(opts, 'beta_2'),   HH = HH + 2 * opts.beta_2 * eye(k); end
    [W, iterW] = doAlIter(HH, VH, W, maxIterW, tolW, 0, maxNumberThreads);
    %fprintf('Learing: %.4E\n\n', norm(V-W*H, 'fro')^2);
    
    HIS.time = [HIS.time toc(startTime)];
    
    HIS.error = [HIS.error norm(V - W'*H, 'fro')^2/2.0];
    
    %if opts.verbose, fprintf('%.10E\n', HIS.error(i+1)); end;
    if isfield(opts, 'mu_1'), HIS.error(i+1) = HIS.error(i+1) + sum(sum(W)) * opts.mu_1; end;
    if isfield(opts, 'mu_2'), HIS.error(i+1) = HIS.error(i+1) + norm(W, 'fro') * opts.mu_2; end;
    if isfield(opts, 'beta_1'), HIS.error(i+1) = HIS.error(i+1) + sum(sum(H)) * opts.beta_1; end;
    if isfield(opts, 'beta_2'), HIS.error(i+1) = HIS.error(i+1) + norm(H, 'fro') * opts.beta_2; end;
    %if opts.verbose, fprintf('%.10E\n', HIS.error(i+1)); end;
     
    [gradW, gradH] = getGrad(V, W', H);
    normGradW = normGrad(gradW, W');
    normGradH = normGrad(gradH, H);
    HIS.gradW = [HIS.gradW normGradW];
    HIS.gradH = [HIS.gradH normGradH];
    HIS.grad = [HIS.grad (normGradW+normGradH)];
    HIS.iterW = [HIS.iterW iterW];
    HIS.iterH = [HIS.iterH iterH];
    HIS.iter = [HIS.iter (iterW*n + iterH*m)/(n+m)];
    if opts.verbose,
        fprintf('%d,\t%.4f,\t%.4E,\t%.4E,\t%.4E,\t%.2E\t%f\t%f,\t%f\n', ...
            i, HIS.time(i+1), HIS.error(i+1), HIS.gradW(i+1), HIS.gradH(i+1), HIS.grad(i+1), HIS.iterW(i+1), HIS.iterH(i+1), HIS.iter(i+1));
    end
end
W=W';
end