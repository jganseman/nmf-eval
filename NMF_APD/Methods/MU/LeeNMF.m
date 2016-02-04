function [W, H, HIS] = LeeNMF(V, opts)
%NNMF - non-negative matrix factorization
%
% [W, H] = nnmf(V, r, iterations)
%
% Input:
%   V   - the matrix to factorize
%   r   - number of basis vectors to generate
%   iterations - number of EM iterations to perform
%
% Results:
%   W   - a set of r basis vectors
%   H   - represenations of the columns of V in 
%         the basis given by W
%
%   Author: David Ross
%   Follows: D.D. Lee & H.S. Seung, "Learning the parts of 
%            objets by non-negative matrix factorization"
%            NATURE, vol 401, 21 Oct. 1999.

%---------------------------------------
% check the input
%---------------------------------------
%[M, N] = size(V);
%---------------------------------------
% Initialization
%---------------------------------------

N = size(V,1); % dimensionality of examples (# rows)

W = opts.W0;
H = opts.H0;


%---------------------------------------
% EM
%---------------------------------------
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

for i = 1:opts.maxIter,
    startTime = tic;
    %% E-Step
    V_over_WH = V ./ (W * H + (V==0));
    W = W .* (V_over_WH * H');
    W = W ./ repmat(sum(W,1), [N 1]);
    
    %% M-Step
    V_over_WH = V ./ (W * H + (V==0));
    H = H .* (W' * V_over_WH);
    
    HIS.time = [HIS.time toc(startTime)];
    HIS.error = [HIS.error norm(V - W*H, 'fro')^2/2.0];
    [gradW, gradH] = getGrad(V, W, H);
    normGradW = normGrad(gradW, W);
    normGradH = normGrad(gradH, H);
    HIS.gradW = [HIS.gradW normGradW];
    HIS.gradH = [HIS.gradH normGradH];
    HIS.grad = [HIS.grad (normGradW+normGradH)];
    HIS.iterW = [HIS.iterW 1];
    HIS.iterH = [HIS.iterH 1];
    if opts.verbose,
        fprintf('%d,\t%.4f,\t%.2f,\t%.2f,\t%.2f,\t%.2E\t%d\n', ...
            i, HIS.time(i+1), HIS.error(i+1), HIS.gradW(i+1), HIS.gradH(i+1), HIS.grad(i+1), HIS.iterW(i+1)+HIS.iterH(i+1));
    end
end
HIS.iters = HIS.iterW + HIS.iterH;
