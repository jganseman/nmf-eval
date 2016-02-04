function [H, HIS] = ActiveSetProject(V, W, opts)
% Nonnegative Matrix Factorization (NMF) via Alolosided Algorithm + Greedy Coordinate Descent
%
% Usage: [H, HIS] = AloProject(V, W, opts)
%
% This software solve the following least squares NMF problem:
%
%	min_{W,H} ||V-WH||_F^2     s.t. W>=0, H>=0
%
% input: 
%		V: the input n by m dense matrix
%       W: the input n by k dense matrix
%		opts.maxIter: maximum number of iterations
%		opts.verbose: 1: compute objective value per iteration.
%                     0: do not compute objective value per iteration. (default)
%
% output: 
%		H: k by m dense matrix
%		HIS.error: objective values.
%		HIS.time : time taken by AloProj.
%       HIS.grad : square of gradient of coordinates that do not satisfy
%       KTT conditions.
%
maxIter = opts.maxIter;

n = size(V, 1);
m = size(V, 2);
k = size(W, 2);

startTime = tic;
HIS.time = toc(startTime);

%[W,gradW,subIterW] = nnlsm(H',A',W',par.nnls_solver);, 
%[H,gradHX,subIterH] = nnlsm(W,A,H,par.nnls_solver);
par.nnls_solver = 'as';
[H, ~, HIS.iter] = nnlsm(W, V, zeros(k, m), par.nnls_solver);

HIS.time = toc(startTime);

HIS.error =  norm(V-W*H, 'fro')^2/2.0;
[~, gradH] = getGrad(V, W, H);
HIS.grad = normGrad(gradH, H);
fprintf('%e, %e, %e, %d\n', HIS.error, HIS.iter, HIS.grad, HIS.time);
end


function [X,grad,iter] = nnlsm(A,B,init,solver)
    switch solver
        case 'bp'
            [X,grad,iter] = nnlsm_blockpivot(A,B,0,init);
        case 'as'
            [X,grad,iter] = nnlsm_activeset(A,B,1,0,init);
    end
end 