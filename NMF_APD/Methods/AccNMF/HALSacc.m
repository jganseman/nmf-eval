% Accelerated hierarchical alternating least squares (HALS) algorithm of
% Cichocki et al. 
%
% See N. Gillis and F. Glineur, "Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization?, 
% Neural Computation 24 (4), pp. 1085-1105, 2012. 
% See http://sites.google.com/site/nicolasgillis/ 
%
% [U,V,e,t] = HALSacc(M,U,V,alpha,delta,maxiter,timelimit)
%
% Input.
%   M              : (m x n) matrix to factorize
%   (U,V)          : initial matrices of dimensions (m x r) and (r x n)
%   alpha          : nonnegative parameter of the accelerated method
%                    (alpha=0.5 seems to work well)
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
% Remark. With alpha = 0, it reduces to the original HALS algorithm.  

function [U, V, HIS] = HALSacc(M, opts)

U = opts.W0;
V = opts.H0;
alpha = opts.alpha;
delta = opts.delta;
maxiter = opts.maxIter;
timelimit = opts.timeLimit;

% Initialization
startTime = tic;
etime = cputime; nM = norm(M,'fro')^2; 
[~,n] = size(M); [m,r] = size(U);
a = 0; e = []; t = []; iter = 0; 

if nargin <= 3, alpha = 0.5; end
if nargin <= 4, delta = 0.1; end
%if nargin <= 5, maxiter = 100; end
%if nargin <= 6, timelimit = 60; end


HIS.error = norm(M - U*V, 'fro')^2/2.0;
% Scaling, p. 72 of the thesis
eit1 = cputime; A = M*V'; B = V*V'; eit1 = cputime-eit1; j = 0;
scaling = sum(sum(A.*U))/sum(sum( B.*(U'*U) )); U = U*scaling;

% Main loop
startTime = tic;
HIS.time = toc(startTime);
[gradW, gradH] = getGrad(M, U, V);
normGradW = normGrad(gradW, U);
normGradH = normGrad(gradH, V);
HIS.gradW = normGradW;
HIS.gradH = normGradH;
HIS.iterW = 0;
HIS.iterH = 0;
HIS.grad = (normGradW + normGradH);
fprintf('%.4E\n', HIS.error);

for i=1:maxiter,
    startTime = tic;
    % Update of U
    if j == 1, % Do not recompute A and B at first pass
        % Use actual computational time instead of estimates rhoU
        eit1 = cputime; A = M*V'; B = V*V'; eit1 = cputime-eit1; 
    end
    j = 1; eit2 = cputime; eps = 1; eps0 = 1;
    [U, iterW] = HALSupdt(M',U',B',A',eit1,alpha,delta); U = U';
    % Update of V
    eit1 = cputime; A = (U'*M); B = (U'*U); eit1 = cputime-eit1;
    eit2 = cputime; eps = 1; eps0 = 1; 
    [V, iterH] = HALSupdt(M,V,B,A,eit1,alpha,delta); 
    % Evaluation of the error e at time t
    if nargout >= 3
        cnT = cputime;
        e = [e ( (nM-2*sum(sum(V.*A))+ sum(sum(B.*(V*V')))) /2)]; 
        etime = etime+(cputime-cnT);
    end
    j = 1;
    HIS.time = [HIS.time toc(startTime)];
    HIS.error = [HIS.error ( (nM-2*sum(sum(V.*A))+ sum(sum(B.*(V*V')))) /2)];
    [gradW, gradH] = getGrad(M, U, V);
    normGradW = normGrad(gradW, U);
    normGradH = normGrad(gradH, V);
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

% Update of V <- HALS(M,U,V)
% i.e., optimizing min_{V >= 0} ||M-UV||_F^2 
% with an exact block-coordinate descent scheme
function [V, iter] = HALSupdt(M,V,UtU,UtM,eit1,alpha,delta)
[r,n] = size(V); 
eit2 = cputime; % Use actual computational time instead of estimates rhoU
cnt = 1; % Enter the loop at least once
eps = 1; eps0 = 1; eit3 = 0;
iter = 0;
while cnt == 1 || (cputime-eit2 < (eit1+eit3)*alpha && eps >= (delta)^2*eps0)
    nodelta = 0; if cnt == 1, eit3 = cputime; end
    iter = iter + 1;
        for k = 1 : r
            deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
            V(k,:) = V(k,:) + deltaV;
            nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta; 
        eit3 = cputime-eit3; 
    end
    eps = nodelta; cnt = 0;
end