function [W,H,errs,vout] = nmf_kl_con(V,r,varargin)
% function [W,H,errs,vout] = nmf_kl_con(V,r,varargin)
%
% Implements Convolutive NMF as described in [1]:
%
%      min D(V||W*H) s.t. W>=0, H>=0
% 
% where V = sum_t W(t) shift(H,t) and the shift function moves H's columns
% t positions to the right (introducing columns of 0s on the left as
% necessary).  Negative values of t shift to the left.
%
% Inputs: (all except V and r are optional and passed in in name-valuepairs)
%   V      [mat]  - Input matrix (n x m)
%   r      [num]  - Rank of the decomposition
%   win    [num]  - Width in columns of each W basis matrix [1]
%   niter  [num]  - Max number of iterations to use [100]
%   thresh [num]  - Number between 0 and 1 used to determine convergence;
%                   the algorithm has considered to have converged when:
%                   (err(t-1)-err(t))/(err(1)-err(t)) < thresh
%                   ignored if thesh is empty [[]]
%   norm_w [num]  - Type of normalization to use for W [1]
%                   Can be (see normalize_W for details):
%                      -2 for (patch-level 2-norm)
%                      -1 for (patch-level 1-norm)
%                       0 for no normalization
%                       2 for (column-level 2-norm)
%                       1 for (column-level 1-norm)
%   norm_h [num]  - Type of normalization to use for rows of H [0]
%                   Can be (see normalize_H for details):
%                      -2 for (matrix-level 2-norm)
%                      -1 for (matrix-level 1-norm)
%                       0 for no normalization
%                       2 for (row-level 2-norm)
%                       1 for (row-level 1-norm)
%   verb   [num]  - Verbosity level (0-3, 0 means silent) [1]
%   W0     [mat]  - Initial W values (n x r) [[]]
%                   empty means initialize randomly
%   H0     [mat]  - Initial H values (r x m) [[]]
%                   empty means initialize randomly
%   W      [mat]  - Fixed value of W (n x r) [[]] 
%                   empty means we should update W at each iteration while
%                   passing in a matrix means that W will be fixed
%   H      [mat]  - Fixed value of H (r x m) [[]] 
%                   empty means we should update H at each iteration while
%                   passing in a matrix means that H will be fixed
%   myeps  [num]  - Small value to add to denominator of updates [1e-20]
%
% Outputs:
%   W      [mat]  - Basis matrix (n x r x win)
%   H      [mat]  - Weight matrix (r x m)
%   errs   [vec]  - Error of each iteration of the algorithm
%
% [1] Smaragdis, P. Non-negative Matrix Factor Deconvolution; 
%     Extraction of Multiple Sound Sources from Monophonic Inputs,
%     International Symposium on ICA and BSS, 2004.
%
% 2010-01-14 Graham Grindlay (grindlay@ee.columbia.edu)

% Copyright (C) 2008-2028 Graham Grindlay (grindlay@ee.columbia.edu)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% do some sanity checks
if min(min(V)) < 0
    error('Matrix entries can not be negative');
end
if min(sum(V,2)) == 0
    error('Not all entries in a row can be zero');
end

[n,m] = size(V);

% process arguments
[win, niter, thresh, norm_w, norm_h, verb, myeps, W0, H0, W, H] = ...
    parse_opt(varargin, 'win', 1, 'niter', 100, 'thresh', [], ...
                        'norm_w', 1, 'norm_h', 0, 'verb', 1, ...
                        'myeps', 1e-20, 'W0', [], 'H0', [], ...
                        'W', [], 'H', []);

% initialize W based on what we got passed
if isempty(W)
    if isempty(W0) % initialize randomly
        W = rand(n,r,win);
    else % use the initial W values we got passed
        W = W0;
    end
    update_W = true;
else % we aren't updating W
    update_W = false;
end

% initialize H based on what we got passed
if isempty(H)
    if isempty(H0) % initialize randomly
        H = rand(r,m);
    else % use the initial H values we got passed
        H = H0;
    end
    update_H = true;
else % we aren't updating H
    update_H = false;
end
                    
if norm_w ~= 0
    % normalize W
    W = normalize_W(W,norm_w);
end

if norm_h ~= 0
    % normalize H
    H = normalize_H(H,norm_h);
end

% preallocate matrix of ones
Onm = ones(n,m);

errs = zeros(niter,1);
for t = 1:niter
    % update W if requested
    if update_W
        for k = 0:win-1
            % approximate V
            R = rec_cnmf(W,H,myeps);
            
            % update W
            W(:,:,k+1) = W(:,:,k+1) .* (((V./R)*shift(H,k)') ./ ...
                                        max(Onm*shift(H,k)', myeps));
        end
        
        if norm_w ~= 0
            % normalize W
            W = normalize_W(W,norm_w);
        end
    end
    
    % update reconstruction
    R = rec_cnmf(W,H,myeps);
   
    % update H if requested
    if update_H
        h = H;
        H = zeros(r,m);
        for k = 0:win-1
            H = H + h .* (W(:,:,k+1)'*shift((V./R),-k)) ./ ...
                         max(W(:,:,k+1)'*Onm, myeps);
        end
        H = H ./ win;
    
        if norm_h ~= 0
            H = normalize_H(H,norm_h);
        end
    end
    
    % update reconstruction
    R = rec_cnmf(W,H,myeps);
    
    % compute I-divergence
    errs(t) = sum(V(:).*log(V(:)./R(:)) - V(:) + R(:));
    
    % display error if asked
    if verb >= 3
        disp(['nmf_kl_con: iter=' num2str(t) ', err=' num2str(errs(t))]);
    end
    
    % check for convergence if asked
    if ~isempty(thresh)
        if t > 2
            if (errs(t-1)-errs(t))/(errs(1)-errs(t-1)) < thresh
                break;
            end
        end
    end
end

% display error if asked
if verb >= 2
    disp(['nmf_kl_con: final_err=' num2str(errs(t))]);
end

% if we broke early, get rid of extra 0s in the errs vector
errs = errs(1:t);

% needed to conform to function signature required by nmf_alg
vout = {};

end

function R = rec_cnmf(W,H,myeps)
% function R = rec_cnmf(W,H,myeps)
%
% Reconstruct a matrix R using Convolutive NMF using W and H matrices.
%

[n, r, win] = size(W);
m = size(H,2);

R = zeros(n,m);
for t = 0:win-1
    R = R + W(:,:,t+1)*shift(H,t);
end

R = max(R,myeps);

end

function O = shift(I, t)
% function O = shift(I, t)
%
% Shifts the columns of an input matrix I by t positions.  
% Zeros are shifted in to new spots.
%

if t < 0
    O = [I(:,-t+1:end) zeros(size(I,1),-t) ];
else
    O = [zeros(size(I,1),t) I(:,1:end-t) ];
end

end
