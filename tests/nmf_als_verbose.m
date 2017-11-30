function [W,H,errs,iterW,iterH,logspectdist]=nmf_als_verbose(X,K,varargin)
% function [W,H,errs,iterW,iterH,logspectdist] = nmf_als_verbose(X,K,Nitsmax,speak)
% co-author: Joachim Ganseman
%
% This function is a copy of [W,H]=nmf_als(X,K,Nitsmax,speak)
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
% ----------
%
% Truncated linear solution to LS -NMF
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% maxiter : Maximum number of iterations to run
% speak   : prints iteration count and changes in connectivity matrix 
%           elements unless speak is 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Lars Kai Hansen, IMM-DTU (c) October 2006

% JGA: Additional arguments to start with a given matrix and do input checks
if min(min(X)) < 0
    error('Matrix entries can not be negative');
end
if min(sum(X,2)) == 0
    error('Not all entries in a row can be zero');
end

% JGA: make arguments into parsable list 
[Nitsmax, speak, thresh, W0, H0] = ...
    parse_opt(varargin, 'Nitsmax', 100, 'speak', 0, 'thresh', 10^(-5), ...
                        'W0', [], 'H0', []);

[D,N]=size(X);
Xscale=sum(sum(X));                    

% JGA: allow for initialization with external matrices
if isempty(W0)
    W = rand(D,K);
else
    W = W0;
end
if isempty(H0)
    H = rand(K,N);
else
    H = H0;
end

% JGA: init additional output variables
errs = zeros(Nitsmax,1);
iterW = zeros(D, K, Nitsmax);     % JGA : addition
iterH = zeros(K, N, Nitsmax);     % JGA : addition
logspectdist = zeros(Nitsmax,1);      % JGA : addition                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print_iter = 50; % iterations between print on screen

%INIT
Rscale=sum(sum(W*H));
sqrnorm=sqrt(Rscale/Xscale);
H=H/sqrnorm;
W=W/sqrnorm;

Xr_old = W*H;

%ITERATE
for n=1:Nitsmax,
    %    W=X*(pinv(H*H')*H)'; % old updates
    W = ((pinv(H*H')*H)*X')';
    W=(W>0).*W+eps;       % JGA: added eps to avoid 0 values (and subsequent NaN)
    W=W./(repmat(sum(W),D,1)+eps); % normalize columns to unit length

    H=(W*pinv(W'*W))'*X;
    H=H.*(H>0)+eps;       % JGA: added eps to avoid 0 values (and subsequent NaN)

    % print to screen
    if (rem(n,print_iter)==0) & speak,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist  = nmf_euclidean_dist(X,W*H);
        errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(n),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < tresh, break, end
    end
    
    %JGA: save intermediate results to output parameters
    Xr = W*H;
    errs(n) = sum(sum((X-Xr).^2)); %EUCLIDEAN!
    iterW(:,:,n) = W;
    iterH(:,:,n) = H;
    % JGA : compute mean log spectral distance
    % first: divide, log square elementwise
    % then: reshape back into a matrix size of the input
    % finally: sum column wise and take square roots, calculate mean value
    tempdist = reshape( (10*log10((X(:)+eps)./(Xr(:)+eps))).^2 , D, N);
    %TODO uncomment this to compute those values after performance testing is finished
    logspectdist(n) = mean(sqrt(sum(tempdist,1)));
    
end
