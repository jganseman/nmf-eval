function [W,H,out_errs,out_iterW,out_iterH,logspectdist] = nmf_cjlin_verbose(V,varargin)
% function [W,H,out_errs,out_iterW,out_iterH,logspectdist] = nmf_cjlin_verbose(X,K,Nitsmax,speak)
% co-author: Joachim Ganseman
%
% This function is a copy of [W,H]=nmf_cjlin(V,Winit,Hinit,tol,timelimit,maxiter)
% from Lars Kai Hansen's NMF_DTU library, with some additional output encoded
% within. The original function header is repeated below.
%
% Additionally, we output the mean log spectral distance between subsequent
% iterations of the algorithm, in order to have a symmetric distance metric
% to evaluate convergence behaviour for all values of beta.
%
% Additionally, we also output the estimated spectrograms between
% subsequent iterations. This allows to run the ISTFT later, and e.g.
% compute signal-based metrics. Access nth estimate trough iterspect(:,:,n).
%
% The parts that were added to the original routine are commented with JGA.


% Copyright (c) 2005-2006 Chih-Jen Lin
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
% 1. Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% 3. Neither name of copyright holders nor the names of its contributors
% may be used to endorse or promote products derived from this software
% without specific prior written permission.
% 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
% NMF by alternative non-negative least squares using projected gradients
% Author: Chih-Jen Lin, National Taiwan University
%
% W,H: output solution
% Winit,Hinit: initial solution
% tol: tolerance for a relative stopping condition
% timelimit, maxiter: limit of time and iterations


% JGA: Additional arguments to start with a given matrix and do input checks
if min(min(V)) < 0
    error('Matrix entries can not be negative');
end
if min(sum(V,2)) == 0
    error('Not all entries in a row can be zero');
end

% JGA: make arguments into parsable list 
[Winit,Hinit,tol,timelimit,maxiter] = ...
    parse_opt(varargin, 'Winit', [], 'Hinit', [], ...
                'tol', eps, 'timelimit', 60, 'maxiter', 100);
%JGA: muted excessive verbose output
%JGA: set tolerance to lowest possible, set timelimit to 1 min

% JGA: init additional output variables
[D,N] = size(V);
K = size(Winit, 2);
out_errs = zeros(maxiter,1);
out_iterW = zeros(D, K, maxiter);     % JGA : addition
out_iterH = zeros(K, N, maxiter);     % JGA : addition
logspectdist = zeros(maxiter,1);      % JGA : addition   




W = Winit; H = Hinit; initt = cputime;

gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;
initgrad = norm([gradW; gradH'],'fro');
%fprintf('Init gradient norm %f\n', initgrad); 
tolW = max(0.001,tol)*initgrad; tolH = tolW;

for iter=1:maxiter
  % stopping condition
  projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
  if projnorm < tol*initgrad || cputime-initt > timelimit
    break;
  end
  
  [W,gradW,iterW] = nlssubprob(V',H',W',tolW,1000); W = W'; gradW = gradW';
  if iterW==1
    tolW = 0.1 * tolW;
  end

  [H,gradH,iterH] = nlssubprob(V,W,H,tolH,1000);
  if iterH==1
    tolH = 0.1 * tolH; 
  end
  
  % JGA: added eps to avoid 0 values (and subsequent NaN)
  W=W+eps;
  H=H+eps;

  %JGA: save intermediate results to output parameters
    Xr = W*H;
    out_errs(iter) = sum(sum((V-Xr).^2)); %EUCLIDEAN!
    out_iterW(:,:,iter) = W;
    out_iterH(:,:,iter) = H;
    % JGA : compute mean log spectral distance
    % first: divide, log square elementwise
    % then: reshape back into a matrix size of the input
    % finally: sum column wise and take square roots, calculate mean value
    tempdist = reshape( (10*log10((V(:)+eps)./(Xr(:)+eps))).^2 , D, N);
    %TODO uncomment this to compute those values after performance testing is finished
    logspectdist(iter) = mean(sqrt(sum(tempdist,1)));
  
  if (iterW==1 && iterH==1 && tolH + tolW < tol*initgrad)
    fprintf('Failed to move\n'); break;
  end
  %if rem(iter,10)==0, fprintf('.'); end
  
end
%fprintf('\nIter = %d Final proj-grad norm %f\n', iter, projnorm);



function [H,grad,iter] = nlssubprob(V,W,Hinit,tol,maxiter)
% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations

H = Hinit; 
WtV = W'*V;
WtW = W'*W; 

alpha = 1; beta = 0.1;
for iter=1:maxiter  
  grad = WtW*H - WtV;
  projgrad = norm(grad(grad < 0 | H >0));
  if projgrad < tol
    break
  end

  % search step size 
  Hn = max(H - alpha*grad, 0); d = Hn-H;
  gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
  if gradd + 0.5*dQd > 0.01*gradd 
    % decrease alpha
    while 1
      alpha = alpha*beta;
      Hn = max(H - alpha*grad, 0); d = Hn-H;
      gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
      if gradd + 0.5*dQd <= 0.01*gradd || alpha < 1e-20    %JGA: fix boolean operators  
        H = Hn; break;
      end
    end 
  else 
    % increase alpha
    while 1
      Hp = Hn;
      alpha = alpha/beta;
      Hn = max(H - alpha*grad, 0); d = Hn-H;
      gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
      if gradd + 0.5*dQd > 0.01*gradd || isequal(Hn,Hp) || alpha > 1e10  %JGA: fix boolean operators    
        H = Hp; alpha = alpha*beta; break;
      end
    end 
  end
end

%if iter==maxiter,
%  fprintf('Max iter in nlssubprob\n');
%end

