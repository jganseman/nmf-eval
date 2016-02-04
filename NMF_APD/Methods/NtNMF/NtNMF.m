function [B, C, HIS] = NtNMF(A, opts, varargin)
% https://www.cs.utexas.edu/~dmkim/Source/software/nnma/nnma.html
% FNMAE : fast nonnegative matrix approximation
% function [B, C, obj] = fnmae(A, maxit, B0, option, varargin)
%
% On output, 
% B and C are nonnegative factors, obj is objective value
%
% On input, 
% A is the target matrix 
% k is rank of B and C.
% maxit is the maximum number of iterations
% B0 is an initial matrix.
% option = {0|1|2}, 0 for Newton, 1 for Quasi-Newton, 
% 2 for CG. Default is 1.
maxit = opts.maxIter;
B0 = opts.W0;
option_ = 0;
global option;

[m,n] = size(A);
option = 1;

if nargin == 4 
  option = option_;
end

% convergence threshold. 
TolFun1 = 1e-2;
TolFun2 = 1e-2;
B = B0;

% initialize
if issparse(B)
  R = qr(B);
else
  R = triu(qr(B));
end
C = R \ (R' \ (B' * A));
C = C + R \ (R' \ (B' * (A - B * C)));
C(C < 0) = 0;

if issparse(C)
  R = qr(C');
else
  R = triu(qr(C'));
end
B = (R \ (R' \ (C * A')))';
B = B + (R \ (R' \ (C * (A' - C' * B'))))';
B(B < 0) = 0;

numCCol = size(C, 2);
numBRow = size(B, 1);
prCol = floor(numCCol / 4);
prRow = floor(numBRow / 4);

% log information
startTime = tic;
HIS.time = toc(startTime);
HIS.error = norm(A - B*C)^2;
[gradW, gradH] = getGrad(A, B, C);
normGradW = normGrad(gradW, B);
normGradH = normGrad(gradH, C);
HIS.gradW = normGradW;
HIS.gradH = normGradH;
HIS.iterW = 0;
HIS.iterH = 0;
HIS.grad = (normGradW + normGradH);

for iter = 1 : maxit,
  startTime = tic;
  if option == 1 || option == 2
	for i = 1 : numCCol,
	  [C(:, i), cObj] = subproc(B, A(:, i), C(:, i), TolFun1);
	  
	  if mod(i, prCol) == 0
    	%fprintf('.');
	  end	
	end
	
	for i = 1 : numBRow,
	  [D(:, i), cObj] = subproc(C', A(i, :)', B(i, :)', TolFun2);
	  
	  if mod(i, prRow) == 0,
          %fprintf('*');
	  end	
	end
  else
	if issparse(B)
	  R = qr(B);
	else
	  R = triu(qr(B));
	end
	for i = 1 : numCCol,
	  [C(:, i), cObj] = subproc(B, A(:, i), C(:, i), TolFun1, R);
	  
	  if mod(i, prCol) == 0,
    	%fprintf('.');
	  end	
	end
	
	if issparse(C)
	  R = qr(C');
	else
	  R = triu(qr(C'));
	end
	for i = 1 : numBRow,
	  [D(:, i), cObj] = subproc(C', A(i, :)', B(i, :)', TolFun2, R);
	  
	  if mod(i, prRow) == 0,
    	%fprintf('*');
	  end	
	end
  end
  
  B = D';
  
  %fprintf('\n');

  if iter < 7
	TolFun1 = TolFun1 * 1e-2;
	TolFun2 = TolFun2 * 1e-2;
  end
  % log informaiton
  HIS.time = [HIS.time toc(startTime)];
  HIS.error = [HIS.error (norm(A - B*C)^2)];
  [gradW, gradH] = getGrad(A, B, C);
  normGradW = normGrad(gradW, B);
  normGradH = normGrad(gradH, C);
  HIS.gradW = [HIS.gradW normGradW];
  HIS.gradH = [HIS.gradH normGradH];
  HIS.grad = [HIS.grad (normGradW + normGradH)];
  iterW = 1;
  iterH = 1;
  HIS.iterW = [HIS.iterW 1];
  HIS.iterH = [HIS.iterH 1];
  %size(HIS.gradW)
  i = iter;
  %size(HIS.time), size(HIS.error), size(HIS.gradW)
  if opts.verbose,
      fprintf('%d,\t%.4f,\t%.2f,\t%.2f,\t%.2f,\t%.2E,\t%d\n', ...
          i, HIS.time(i+1), HIS.error(i+1), HIS.gradW(i+1), HIS.gradH(i+1), HIS.grad(i+1), (iterW+iterH));
  end
end

obj = 0.5 * norm(A - B * C, 'fro')^2;

function [x, f] = subproc(A, b, x, TolFun, R, varargin) 
  global option Aglb bglb invHess Rglb AtAu nAu deltaX;

  % initialize parameters
  Aglb = A;
  bglb = b;
  [numRows, numCols] = size(Aglb);

  f = func(x);                                % objective value
  grad = gradf(x);                            % gradient
  done = false;
  iter = 0;
  updateHess = 1;

  if nargin == 4 && option ~= 2
    updateHess = 1;
	invHess = eye(numCols);
	findSrch = @(grad)(findSrchQN(grad));
  elseif option == 2
	updateHess = 2;
	nAu = 1;
	AtAu = zeros(numCols, 1);
	deltaX = zeros(numCols, 1);
	findSrch = @(grad)(findSrchCG(grad));
  elseif option == 0
	updateHess = 0;
	Rglb = R;
	findSrch = @(grad)(findSrchN(grad));
  end
  
  % Main loop
  while ~done
    iter = iter + 1;
    xOld = x; 
	zx = find(x == 0);
	gp = find(grad > 0);
	gp = intersect(zx, gp);
	grad(gp) = 0;                  % \nabla f(\vec{y}^k)

    srchDir = findSrch(grad);
	ns = find(srchDir < 0);
	ns = intersect(zx, ns);
	srchDir(gp) = 0;               % d
	srchDir(ns) = 0;               % \bar{d}

	A = Aglb;
	A(:,gp) = 0;
	Ad = A * srchDir;
	alphaLnSrch = -Ad' * ((A * x) - b) / (Ad' * Ad);

	if norm(Ad) <= max(TolFun, 1e-8) || alphaLnSrch <= max(TolFun, 1e-8)
	  break;
	end
	
	x = x + alphaLnSrch * srchDir;
	x(x < 1e-8) = 0;
	deltaX = x - xOld;

	grad = gradf(x);
	fNew = func(x);

	if abs(f - fNew) <= max(TolFun, 1e-8)
	  break;
	else
	  f = fNew;
	end
	
    % Updates the quasi-Newton matrix that approximates
    % the inverse to the Hessian.
	if updateHess == 1
	  Au = A * deltaX;
	  nAu = norm(Au)^2;
	  AtAu = A' * Au;
	  DAtAu = invHess * AtAu;
	  invHess = invHess + ((1 + AtAu' * DAtAu / nAu) * ...
						   deltaX * deltaX' - ...
						   (DAtAu * deltaX' + deltaX * DAtAu')) / nAu;
	elseif updateHess == 2
	  Au = A * deltaX;
	  nAu = norm(Au)^2;
	  AtAu = A' * Au;
	end
  end % of while

  f = func(x);                                % objective value
  
%--------------------------------------------------------------------------
% Objective functions and gradients
%
function fc = func(x)
  global Aglb bglb;

  Ax = Aglb * x;
  fc = norm(Ax - bglb)^2;
  
function gc = gradf(x)
  global Aglb bglb;

  Ax = Aglb * x;
  gc = (Aglb' * (Ax - bglb));
  
%--------------------------------------------------------------------------
% Find search direction
%
function srchDir = findSrchQN(grad)
  global invHess;
  
  srchDir = -invHess * grad;

function srchDir = findSrchN(grad)
  global Rglb;

  srchDir = Rglb \ (Rglb' \ -grad);
  
function srchDir = findSrchCG(grad)
  global AtAu nAu deltaX;

  srchDir = -grad + ((grad' * AtAu) / nAu) * deltaX;
