% Accelerated Multiplicative Updates of Lee and Seung
%
% See N. Gillis and F. Glineur, "Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization?, 
% Neural Computation 24 (4), pp. 1085-1105, 2012. 
% See http://sites.google.com/site/nicolasgillis/ 
%
% [U,V,e,t] = MUacc(M,U,V,alpha,delta,maxiter,timelimit)
%
% Input.
%   M              : (m x n) matrix to factorize
%   (U,V)          : initial matrices of dimensions (m x r) and (r x n)
%   alpha          : nonnegative parameter of the accelerated method
%                    (alpha=2 seems to work well)
%   delta          : parameter to stop inner iterations when they become
%                    inneffective (delta=0.1 seems to work well). 
%   maxiter        : maximum number of iterations
%   timelimit      : maximum time alloted to the algorithm
%
% Output.
%   (U,V)    : nonnegative matrices s.t. UV approximate M
%   (e,t)    : error and time after each iteration, 
%               can be displayed with plot(t,e)
%
% Remark. With alpha = 0, it reduces to the original MU of Lee and Seung. 

function [U,V,e,t] = MUacc(M,U,V,alpha,delta,maxiter,timelimit)

% Initialization
startTime = tic;
etime = cputime; nM = norm(M,'fro')^2; 
[m,n] = size(M); [m,r] = size(U);
if issparse(M), K = sum(M(:) ~= 0); else K = m*n; end
rhoU = 1+(K+n*r)/(m*(r+1)); 
rhoV = 1+(K+m*r)/(n*(r+1)); 
a = 0; e = []; t = [];  iter = 0; 

if nargin <= 3, alpha = 2; end
if nargin <= 4, delta = 0.1; end
if nargin <= 5, maxiter = 100; end
if nargin <= 6, timelimit = 60; end

% Main loop
while iter < maxiter && cputime-etime <= timelimit 
    % Update of U
    A = (M*V'); B = (V*V'); 
    eps = 1; eps0 = 1; j = 1;
    while j <= floor(1+rhoU*alpha) && eps >= delta*eps0
        U0 = U; 
        U = max(1e-16,U.*(A./(U*B))); 
        if j == 1
            eps0 = norm(U0-U,'fro'); 
        end
        eps = norm(U0-U,'fro');  j = j+1;
    end
    % Update of V
    A = (U'*M); B = (U'*U); 
    eps = 1; eps0 = 1; j = 1; 
    while j <= floor(1+rhoV*alpha)  &&  eps >= delta*eps0
        V0 = V;
        V = max(1e-16,V.*(A./(B*V))); 
        if j == 1
            eps0 = norm(V0-V,'fro'); 
        end
        eps = norm(V0-V,'fro');  j = j+1;
    end
    % Evaluation of the error e at time t 
    if nargout >= 3
        cnT = cputime;
        e = [e ((nM-2*sum(sum(V.*A))+ sum(sum(B.*(V*V')))) / 2 )];
        etime = etime+(cputime-cnT);    
        t = [t toc(startTime)]; 
    end
    iter = iter+1; 
    disp(sprintf('%d, %.20f', iter, e(iter)));
end