function op = proj_nuclear( q )

%PROJ_NUCLEAR  Projection onto the set of matrices with nuclear norm less than or equal to q.
%    OP = PROJ_NUCLEAR( q ) returns a function that implements the
%    indicator for matrices with nuclear norm less than q.
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied, 
%    it must be a positive real scalar.
%    This function is like proj_psdUTrace.m but does not assume
%    that inputs are square Hermitian positive semidefinite matrices.
%
% This version uses a dense svd decomposition; future versions
% of TFOCS will take advantage of low-rank and/or sparse structure.
% Dual function: prox_spectral.m
% See also prox_spectral, proj_psdUTrace.m, proj_simplex.m, prox_nuclear.m

% June 26, 2012

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
q = proj_simplex( q );
op = @(varargin)proj_nuclear_q( q, varargin{:} );

function [ v, X ] = proj_nuclear_q( eproj, X, t )

VECTORIZE = false;
% Input must be a matrix

% In proj_psdUTrace.m, we allow vector inputs because
%   we can reshape them (since we have square matrices).
%   For nuclear norm, inputs can be rectangular matrices,
%   so we do not know how to correctly reshape a vector in
%   this case.

% if size(X,1) ~= size(X,2)
%     n = sqrt(length(X));
%     X = reshape(X, n, n );
%     VECTORIZE   = true;
% end
v = 0;
if nargin > 2 && t > 0,
    [U,D,V] = svd(X,'econ');
    [dum,D] = eproj(diag(D),t);
    tt = D > 0;
    X  = U(:,tt)*diag(D(tt))*V(:,tt)';
%     if VECTORIZE, X = X(:); end
elseif any(svd(X)<0),
    v = Inf;
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.

