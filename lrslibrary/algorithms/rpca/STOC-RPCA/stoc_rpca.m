% Stochastic optimization for the robust PCA
% Input:
%       D: [m x n] data matrix, m - ambient dimension, n - samples number
%       lambda1, lambda2: trade-off parameters
%       nrank: the groundtruth rank of the data
% Output:
%       L: [m x r] the basis of the subspace
%       R: [r x n] the coefficient of the samples on the basis L
%       E: sparse error
%
% copyright (c) Jiashi Feng (jshfeng@gmail.com)
%
%{
lambda1 = 1/sqrt(max(size(M)));
lambda2 = lambda1;
nrank = size(M,2);
[L,S] = stoc_rpca(M, lambda1, lambda2, nrank);
show_2dvideo(M,m,n);
show_2dvideo(L,m,n);
show_2dvideo(S,m,n);
%}
% [L,R,E] = stoc_rpca(D, lambda1, lambda2, nrank)
function [Lr,S] = stoc_rpca(D, lambda1, lambda2, nrank)
  %% initialization
  [ndim, nsample] = size(D);
  %L = cell(nsample+1,1);
  %L{1} = rand(ndim,nrank);
  L = rand(ndim,nrank);
  A = zeros(nrank,nrank);
  B = zeros(ndim,nrank);
  %R = zeros(nrank,nsample);
  S = zeros(ndim,nsample);
  Lr = zeros(ndim,nsample);
  %% online optimization
  for t = 1:nsample
    disp(t);
    %   tic;
    z = D(:,t);
    %   tic;
    %[r,s] = solve_proj2(z,L{t},lambda1,lambda2);
    [r,s] = solve_proj2(z,L,lambda1,lambda2);
    %   tused = toc;
    %   fprintf('elapsed time for projection %f secs\n',tused);
    %   R(:,t) = r;
    S(:,t) = s;
    A = A + r*r';
    B = B + (z-s)*r';
    %   L{t+1} = update_col(L{t},A,B,lambda1/nsample);
    L = update_col(L,A,B,lambda1/nsample);
    Lr(:,t) = L * r;
    %   disp(t);
    %   toc;
  end
end

function L = update_col(L,A,B,lambda1)
  [junk,ncol] = size(L);
  A = A + lambda1*eye(ncol,ncol);
  for j = 1:ncol
    bj = B(:,j);
    lj = L(:,j);
    aj = A(:,j);
    temp = (bj-L*aj)/A(j,j) + lj;
    L(:,j) = temp/max(norm(temp),1);
  end
end