function [W, H, HIS] = NMF_GCD(V, opts)
% Nonnegative Matrix Factorization (NMF) via Greedy Coordinate Descent
%
% Usage: [W H objGCD timeGCD] = NMF_GCD(V, k, maxiter, Winit, Hinit, trace)
%
% This software solve the following least squares NMF problem:
%
%	min_{W,H} ||V-WH||_F^2     s.t. W>=0, H>=0
%
% input: 
%		V: the input n by m dense matrix
%		k: the specified rank
%		maxiter: maximum number of iterations
%		Winit: initial of W (n by k dense matrix)
%		Hinit: initial of H (k by m dense matrix)
%		trace: 1: compute objective value per iteration.
%			   0: do not compute objective value per iteration. (default)
%
% output: 
%		NMF_GCD will output nonnegative matrices W, H, such that WH is an approximation of V
%		W: n by k dense matrix
%		H: k by m dense matrix
%		objGCD: objective values. 
%		timeGCD: time taken by GCD. 
%
maxiter = opts.maxIter;
Winit = opts.W0;
Hinit = opts.H0;
tol = opts.tolerance;

n = size(V,1);
m = size(V,2);
W = Winit;
H = Hinit;
k = size(Winit, 2);
%% Stopping tolerance for subproblems


total = 0;

startTime = tic;
HIS.time = toc(startTime);
HIS.error = norm(V - W*H, 'fro')^2/2.0;
[gradW, gradH] = getGrad(V, W, H);
normGradW = normGrad(gradW, W);
normGradH = normGrad(gradH, H);
HIS.gradW = normGradW;
HIS.gradH = normGradH;
HIS.iterW = 0;
HIS.iterH = 0;
HIS.grad = (normGradW+normGradH);

for iter = 1:maxiter,
    startTime = tic;

	% update variables of H
	WV = W'*V;
	WW = W'*W;
	GH = -(WV-WW*H);
	
	[Hnew, iterH] = doiter(GH, WW,H,tol, k^2); % Coordinate descent updates for H
	H = Hnew;

	% update variables of W
	VH = V*H';
	HH = H*H'; % Hessian of each row of W
	
	GW = -(VH - W*HH); % gradient of W

	[Wnew, iterW] = doiter(GW', HH, W',tol, k^2); % Coordinate descent updates for W%
	W = Wnew';
    
    i = iter;
    
    HIS.time = [HIS.time toc(startTime)];
    HIS.error = [HIS.error norm(V - W*H, 'fro')^2/2.0];
    [gradW, gradH] = getGrad(V, W, H);
    normGradW = normGrad(gradW, W);
    normGradH = normGrad(gradH, H);
    HIS.gradW = [HIS.gradW normGradW];
    HIS.gradH = [HIS.gradH normGradH];
    HIS.grad = [HIS.grad (normGradW+normGradH)];
    HIS.iterW = [HIS.iterW iterW];
    HIS.iterH = [HIS.iterH iterH];
    if opts.verbose,
        fprintf('%d,\t%.4f,\t%.2f,\t%.2f,\t%.2f,\t%.2E\t%d,\t%d\n', ...
            i, HIS.time(i+1), HIS.error(i+1), HIS.gradW(i+1), HIS.gradH(i+1), HIS.grad(i+1), HIS.iterW(i+1), HIS.iterH(i+1));
    end
end
