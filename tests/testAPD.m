% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% may 2017

% This demo tests different implementations of NMF, as found in the 
% NMF_APD library, with audio data.

% For the STFT, we choose the corrected version from CATbox, which is available
% in github.com/jganseman/nsgt_eval . This is the fastest STFT among the 4 
% we tested, and also the one with the smallest reconstruction error when
% window-corrected (option 'smooth').

% Phase is ignored in the NMF; the inverse STFT uses the original phase data.

% For separation instead of signal reconstruction, BSS_EVAL 3.0 is used.

%%
clear;
addpath(genpath('../'));     % add subdirectories to path. 
    % Make sure you are in the '/tests' folder, not in the nmf-eval root!
    
%% Parameters

fftsize = 1024;
hopsize = fftsize / 4;      % 75% is a must if STFT coefs are going to be changed
bins_per_octave = 48;       % make this a multiple of 12 for music

inputfile = 'pianoclip4notes.wav';
maxiter = 100;
nrcomponents = 4;

myverbose = 1;

disp('--- Reading input file, monofy and zero-pad. ---')

% Read file
datadir = getDataDirectory();       % directory with example files
[origMix,samplerate] = wavread([datadir inputfile]);

% if stereo, make mono
if size(origMix, 2) > 1
 origMix = sum(origMix, 2) ./2 ; 
end

% Because the signal does not start with 0 but our windows do, we zero-pad
% the beginning as to not lose any information.
% TODO: can be done with less padding: See Hodgkinson STFT implementation.
input = [ zeros(fftsize,1) ; origMix ; zeros(fftsize,1) ];
% afterwards, chop as follows:   
% InvOrigMix = InvOrigMix( fftsize +1 : fftsize +length(origMix) );

%%

disp('--- Computing corrected CATbox STFT with periodic Hann window ---');  
% parameters: signal, window, overlap, fftsize

MySTFT = stft_catbox(input, hann(fftsize, 'periodic'), fftsize-hopsize, fftsize);
% parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'

disp('--- Doing inverse STFT with periodic Hann window (smoothed) ---');

InvMySTFT = istft_catbox(MySTFT, fftsize / hopsize, fftsize, 'smooth')';
InvMySTFT = InvMySTFT(fftsize+1:fftsize+length(origMix));

rec_err = norm(origMix-InvMySTFT)/norm(origMix);
fprintf('      Normalized Reconstruction Error: %e \n',rec_err);

%% separate modulus and phase

MySTFTabs = abs(MySTFT)+eps;
phase = MySTFT./MySTFTabs;
[numrows,numcols] = size(MySTFTabs);

%% TEST APD's AccNMF MUacc
% N. Gillis and F. Glineur, "Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization", 
% Neural Computation 24 (4), pp. 1085-1105, 2012.

% [U,V,e,t] = MUacc(M,U,V,alpha,delta,maxiter,timelimit)
%   M              : (m x n) matrix to factorize
%   (U,V)          : initial matrices of dimensions (m x r) and (r x n)
%   alpha          : nonnegative parameter of the accelerated method
%                    (alpha=2 seems to work well)
%   delta          : parameter to stop inner iterations when they become
%                    inneffective (delta=0.1 seems to work well). 
%   maxiter        : maximum number of iterations
%   timelimit      : maximum time alloted to the algorithm
%   (U,V)    : nonnegative matrices s.t. UV approximate M
%   (e,t)    : error and time after each iteration, 
%               can be displayed with plot(t,e)

[m,n] = size(MySTFTabs);
U0 = rand(m,nrcomponents); V0 = rand(nrcomponents,n);
alpha = 2;      % using default values here
delta = 0.1;       % ties into stopping criterion. check function.
timelimit = 10000000;       % make this so high that it barely counts.

tic;
[ W, H, errpath, timepath ] = MUacc ( MySTFTabs, U0, V0, alpha, delta, maxiter, timelimit);
apd_muacc_time=toc;
apd_muacc = W*H;

fprintf('      APD-MUacc Done in time: %f \n',apd_muacc_time);

%% EVALUATE APD's AccNMF MUacc
%inverse transform and cut to size
apd_muacc_newsig = istft_catbox(apd_muacc.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_muacc_newsig = apd_muacc_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_muacc_W_neg = min(min(W))<0;
apd_muacc_H_neg = min(min(H))<0;
fprintf('      APD-MUacc spectra w/ neg values?: %d \n',apd_muacc_W_neg);
fprintf('      APD-MUacc coeffs w/ neg values?: %d \n',apd_muacc_H_neg);

%compute reconstruction error
apd_muacc_rec_err = norm(origMix-apd_muacc_newsig)/norm(origMix);
fprintf('      APD-MUacc Normalized Reconstruction Error: %e \n',apd_muacc_rec_err);

%compute log spect distance
apd_muacc_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_muacc(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-MUacc Final Log Spectral Distance: %e \n',apd_muacc_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_muacc_SDR, apd_muacc_SIR, apd_muacc_SAR] = bss_eval_sources(apd_muacc_newsig',origMix');
fprintf('      APD-MUacc SDR: %f \t SIR: %f \t SAR: %f \n',apd_muacc_SDR, apd_muacc_SIR, apd_muacc_SAR);

disp('--- Finished ---')


%% TEST APD's AccNMF PGLINacc
% function [W,H,HIS] = PGLINacc(V, opts)

% Modified from: Lin, Projected Gradient Methods for Nonnegative Matrix Factorization,
% Neural Computation, 19, p. 2756-2779, 2007, MIT press. 
% W,H: output solution
% Winit,Hinit: initial solution
% alpha: parameter for acceleration
% delta: parameter for stopping inner iterations
% timelimit, maxiter: limit of time and iterations

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.alpha = 0.5;      % using default values. WARNING: differs from MUacc
opts.delta = 0;       % ties into stopping criterion? check function.
opts.timeLimit = 10000000;       % make this so high that it barely counts.
opts.maxIter = maxiter;
opts.verbose = myverbose;

tic;
[ W, H, HIS ] = PGLINacc ( MySTFTabs, opts);
apd_pglinacc_time=toc;
apd_pglinacc = W*H;

fprintf('      APD-PGLINacc Done in time: %f \n',apd_pglinacc_time);

%% EVALUATE APD's AccNMF PGLINacc
%inverse transform and cut to size
apd_pglinacc_newsig = istft_catbox(apd_pglinacc.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_pglinacc_newsig = apd_pglinacc_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_pglinacc_W_neg = min(min(W))<0;
apd_pglinacc_H_neg = min(min(H))<0;
fprintf('      APD-PGLINacc spectra w/ neg values?: %d \n',apd_pglinacc_W_neg);
fprintf('      APD-PGLINacc coeffs w/ neg values?: %d \n',apd_pglinacc_H_neg);

%compute reconstruction error
apd_pglinacc_rec_err = norm(origMix-apd_pglinacc_newsig)/norm(origMix);
fprintf('      APD-PGLINacc Normalized Reconstruction Error: %e \n',apd_pglinacc_rec_err);

%compute log spect distance
apd_pglinacc_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_pglinacc(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-PGLINacc Final Log Spectral Distance: %e \n',apd_pglinacc_logspectdist);


% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_pglinacc_SDR, apd_pglinacc_SIR, apd_pglinacc_SAR] = bss_eval_sources(apd_pglinacc_newsig',origMix');
fprintf('      APD-PGLINacc SDR: %f \t SIR: %f \t SAR: %f \n',apd_pglinacc_SDR, apd_pglinacc_SIR, apd_pglinacc_SAR);

disp('--- Finished ---')

%% TEST APD's AccNMF HALSacc
% function [W,H,HIS] = HALSacc(V, opts)
% (based on HALS by Cichocki)
%   M              : (m x n) matrix to factorize
%   (U,V)          : initial matrices of dimensions (m x r) and (r x n)
%   alpha          : nonnegative parameter of the accelerated method
%                    (alpha=0.5 seems to work well)
%   delta          : parameter to stop inner iterations when they become
%                    inneffective (delta=0.1 seems to work well). 
%   maxiter        : maximum number of iterations
%   timelimit      : maximum time alloted to the algorithm
%   (U,V)    : nonnegative matrices s.t. UV approximate M
%   (e,t)    : error and time after each iteration, 
%               can be displayed with plot(t,e)
% Remark. With alpha = 0, it reduces to the original HALS algorithm. 

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.alpha = 0.5;      % using default values. WARNING: differs from MUacc
opts.delta = 0.1;       % ties into stopping criterion? check function.
opts.timeLimit = 10000000;       % make this so high that it barely counts.
opts.maxIter = maxiter;
opts.verbose = myverbose;

tic;
[ W, H, HIS ] = HALSacc ( MySTFTabs, opts);
apd_halsacc_time=toc;
apd_halsacc = W*H;

fprintf('      APD-HALSacc Done in time: %f \n',apd_halsacc_time);

%% EVALUATE APD's AccNMF HALSacc
%inverse transform and cut to size
apd_halsacc_newsig = istft_catbox(apd_halsacc.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_halsacc_newsig = apd_halsacc_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_halsacc_W_neg = min(min(W))<0;
apd_halsacc_H_neg = min(min(H))<0;
fprintf('      APD-HALSacc spectra w/ neg values?: %d \n',apd_halsacc_W_neg);
fprintf('      APD-HALSacc coeffs w/ neg values?: %d \n',apd_halsacc_H_neg);

%compute reconstruction error
apd_halsacc_rec_err = norm(origMix-apd_halsacc_newsig)/norm(origMix);
fprintf('      APD-HALSacc Normalized Reconstruction Error: %e \n',apd_halsacc_rec_err);

%compute log spect distance
apd_halsacc_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_halsacc(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-HALSacc Final Log Spectral Distance: %e \n',apd_halsacc_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_halsacc_SDR, apd_halsacc_SIR, apd_halsacc_SAR] = bss_eval_sources(apd_halsacc_newsig',origMix');
fprintf('      APD-HALSacc SDR: %f \t SIR: %f \t SAR: %f \n',apd_halsacc_SDR, apd_halsacc_SIR, apd_halsacc_SAR);

disp('--- Finished ---')


%% TEST APD's ActiveSet 
% Jingu Kim (impl) after Hyunsoo Kim and Haesun Park, Nonnegative Matrix Factorization Based on Alternating Nonnegativity Constrained Least Squares and Active Set Method,
% SIAM Journal on Matrix Analysis and Applications, 2008, 30, 713-730

% function [W,H,his] = ActiveSet(A, opts, varargin)
% note: function was probably adapted by NMF_APD author to facilitate comparison. 
% use opts to have initial matrices W0 and H0. For the rest use varargin:
%       TYPE : 'plain' to use formulation (1)
%               'regularized' to use formulation (2)
%               'sparse' to use formulation (3)
%               Default is 'regularized', which is recommended for quick application testing unless 'sparse' or 'plain' is explicitly needed.
%               If sparsity is needed for 'W' factor, then apply this function for the transpose of 'A' with formulation (3).
%                      Then, exchange 'W' and 'H' and obtain the transpose of them.
%               Imposing sparsity for both factors is not recommended and thus not included in this software.
%        NNLS_SOLVER : 'bp' to use the algorithm in [1]
%                      'as' to use the algorithm in [2]
%                      Default is 'bp', which is in general faster.
%        ALPHA : Parameter alpha in the formulation (2) or (3). 
%                Default is the average of all elements in A. No good justfication for this default value, and you might want to try other values.
%        BETA : Parameter beta in the formulation (2) or (3).
%               Default is the average of all elements in A. No good justfication for this default value, and you might want to try other values.
%        MAX_ITER : Maximum number of iterations. Default is 100.
%        MIN_ITER : Minimum number of iterations. Default is 20.
%        MAX_TIME : Maximum amount of time in seconds. Default is 100,000.
%        W_INIT : (m x k) initial value for W.
%        H_INIT : (k x n) initial value for H.
%        TOL : Stopping tolerance. Default is 1e-3. If you want to obtain a more accurate solution, decrease TOL and increase MAX_ITER at the same time.
%        VERBOSE : 0 (default) - No debugging information is collected.
%                  1 (debugging purpose) - History of computation is returned by 'HIS' variable.
%                  2 (debugging purpose) - History of computation is additionally printed on screen.

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);

tic;
[ W, H, HIS ] = ActiveSet(MySTFTabs, opts, 'NNLS_SOLVER', 'as', 'MAX_ITER', maxiter, 'TOL', eps, 'VERBOSE', myverbose*0);
% BUG: fails when verbose > 0!
apd_actset_time=toc;
apd_actset = W*H;

fprintf('      APD-ActiveSet Done in time: %f \n',apd_actset_time);

%% EVALUATE APD's ActiveSet
%inverse transform and cut to size
apd_actset_newsig = istft_catbox(apd_actset.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_actset_newsig = apd_actset_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_actset_W_neg = min(min(W))<0;
apd_actset_H_neg = min(min(H))<0;
fprintf('      APD-ActiveSet spectra w/ neg values?: %d \n',apd_actset_W_neg);
fprintf('      APD-ActiveSet coeffs w/ neg values?: %d \n',apd_actset_H_neg);

%compute reconstruction error
apd_actset_rec_err = norm(origMix-apd_actset_newsig)/norm(origMix);
fprintf('      APD-ActiveSet Normalized Reconstruction Error: %e \n',apd_actset_rec_err);

%compute log spect distance
apd_actset_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_actset(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-ActiveSet Final Log Spectral Distance: %e \n',apd_actset_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_actset_SDR, apd_actset_SIR, apd_actset_SAR] = bss_eval_sources(apd_actset_newsig',origMix');
fprintf('      APD-ActiveSet SDR: %f \t SIR: %f \t SAR: %f \n',apd_actset_SDR, apd_actset_SIR, apd_actset_SAR);

disp('--- Finished ---')


%% TEST APD's BlockPivot 
% Jingu Kim and Haesun Park, Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons,
%                 In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM'08), 353-362, 2008

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);

tic;
[ W, H, HIS ] = Blockpivot(MySTFTabs, opts, 'NNLS_SOLVER', 'bp', 'MAX_ITER', maxiter, 'TOL', eps, 'VERBOSE', myverbose*0);
% BUG: fails when verbose > 0!
apd_blockp_time=toc;
apd_blockp = W*H;

fprintf('      APD-BlockPivot Done in time: %f \n',apd_blockp_time);

%% EVALUATE APD's BlockPivot
%inverse transform and cut to size
apd_blockp_newsig = istft_catbox(apd_blockp.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_blockp_newsig = apd_blockp_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_blockp_W_neg = min(min(W))<0;
apd_blockp_H_neg = min(min(H))<0;
fprintf('      APD-BlockPivot spectra w/ neg values?: %d \n',apd_blockp_W_neg);
fprintf('      APD-BlockPivot coeffs w/ neg values?: %d \n',apd_blockp_H_neg);

%compute reconstruction error
apd_blockp_rec_err = norm(origMix-apd_blockp_newsig)/norm(origMix);
fprintf('      APD-BlockPivot Normalized Reconstruction Error: %e \n',apd_blockp_rec_err);

%compute log spect distance
apd_blockp_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_blockp(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-BlockPivot Final Log Spectral Distance: %e \n',apd_blockp_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_blockp_SDR, apd_blockp_SIR, apd_blockp_SAR] = bss_eval_sources(apd_blockp_newsig',origMix');
fprintf('      APD-BlockPivot SDR: %f \t SIR: %f \t SAR: %f \n',apd_blockp_SDR, apd_blockp_SIR, apd_blockp_SAR);

disp('--- Finished ---')


%% TEST APD's AloExactNMF 
% Hsieh and Dillon? Possibly adapted by Khuong.
% Nonnegative Matrix Factorization (NMF) via Anti-lopsided + Greedy Coordinate Descent
% NOTE: the code from NMF_APD library is definitely not the original Hsieh/Dillon code
% function [W, H, HIS] = AloExactNMF(V, opts)
%		V: the input n by m dense matrix
%		k: the specified rank
%		opts.maxIter: maximum number of iterations
%		opts.W0: initial of W (n by k dense matrix)
%		opts.H0: initial of H (k by m dense matrix)
%		opts.verbose: 1: compute objective value per iteration.
%			   0: do not compute objective value per iteration. (default)
%		NMF_GCD will output nonnegative matrices W, H, such that WH is an approximation of V
%		W: n by k dense matrix
%		H: k by m dense matrix
%		HIS.errors: objective values verus iterations. 
%		HIS.time: time verus iterations. 

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.verbose = myverbose; 
opts.maxIter = maxiter;
opts.maxNumberThreads = 1;

tic;
[ W, H, HIS ] = AloExactNMF(MySTFTabs, opts);
% BUG: fails when verbose > 0!
apd_aloexact_time=toc;
apd_aloexact = W*H;

fprintf('      APD-AloExact Done in time: %f \n',apd_aloexact_time);

%% EVALUATE APD's AloExactNMF
%inverse transform and cut to size
apd_aloexact_newsig = istft_catbox(apd_aloexact.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_aloexact_newsig = apd_aloexact_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_aloexact_W_neg = min(min(W))<0;
apd_aloexact_H_neg = min(min(H))<0;
fprintf('      APD-AloExact spectra w/ neg values?: %d \n',apd_aloexact_W_neg);
fprintf('      APD-AloExact coeffs w/ neg values?: %d \n',apd_aloexact_H_neg);

%compute reconstruction error
apd_aloexact_rec_err = norm(origMix-apd_aloexact_newsig)/norm(origMix);
fprintf('      APD-AloExact Normalized Reconstruction Error: %e \n',apd_aloexact_rec_err);

%compute log spect distance
apd_aloexact_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_aloexact(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-AloExact Final Log Spectral Distance: %e \n',apd_aloexact_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_aloexact_SDR, apd_aloexact_SIR, apd_aloexact_SAR] = bss_eval_sources(apd_aloexact_newsig',origMix');
fprintf('      APD-AloExact SDR: %f \t SIR: %f \t SAR: %f \n',apd_aloexact_SDR, apd_aloexact_SIR, apd_aloexact_SAR);

disp('--- Finished ---')


%% TEST APD's FCDM NMF 
% Hsieh and Dillon. Possibly adapted by Khuong.
% Nonnegative Matrix Factorization (NMF) via Anti-lopsided + Greedy Coordinate Descent
% function [W, H, HIS] = NMF_GCD(V, opts)
%		V: the input n by m dense matrix
%		k: the specified rank
%		maxiter: maximum number of iterations
%		Winit: initial of W (n by k dense matrix)
%		Hinit: initial of H (k by m dense matrix)
%		trace: 1: compute objective value per iteration.
%			   0: do not compute objective value per iteration. (default)
%		NMF_GCD will output nonnegative matrices W, H, such that WH is an approximation of V
%		W: n by k dense matrix
%		H: k by m dense matrix
%		objGCD: objective values. 
%		timeGCD: time taken by GCD. 

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.verbose = myverbose; 
opts.maxIter = maxiter;
opts.tolerance = eps;

tic;
[ W, H, HIS ] = NMF_GCD(MySTFTabs, opts);
apd_fcdm_time=toc;
apd_fcdm = W*H;

fprintf('      APD-FCDM Done in time: %f \n',apd_fcdm_time);

%% EVALUATE APD's FCDM NMF
%inverse transform and cut to size
apd_fcdm_newsig = istft_catbox(apd_fcdm.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_fcdm_newsig = apd_fcdm_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_fcdm_W_neg = min(min(W))<0;
apd_fcdm_H_neg = min(min(H))<0;
fprintf('      APD-FCDM spectra w/ neg values?: %d \n',apd_fcdm_W_neg);
fprintf('      APD-FCDM coeffs w/ neg values?: %d \n',apd_fcdm_H_neg);

%compute reconstruction error
apd_fcdm_rec_err = norm(origMix-apd_fcdm_newsig)/norm(origMix);
fprintf('      APD-FCDM Normalized Reconstruction Error: %e \n',apd_fcdm_rec_err);

%compute log spect distance
apd_fcdm_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_fcdm(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-FCDM Final Log Spectral Distance: %e \n',apd_fcdm_logspectdist);


% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_fcdm_SDR, apd_fcdm_SIR, apd_fcdm_SAR] = bss_eval_sources(apd_fcdm_newsig',origMix');
fprintf('      APD-FCDM SDR: %f \t SIR: %f \t SAR: %f \n',apd_fcdm_SDR, apd_fcdm_SIR, apd_fcdm_SAR);

disp('--- Finished ---')


%% TEST APD's MU (Ross) 
% Note: code was likely adapted by Khuong!
% function [W, H, HIS] = LeeNMF(V, opts)
%   V   - the matrix to factorize
%   r   - number of basis vectors to generate
%   iterations - number of EM iterations to perform
%   W   - a set of r basis vectors
%   H   - represenations of the columns of V in the basis given by W

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.verbose = myverbose; 
opts.maxIter = maxiter;

tic;
[ W, H, HIS ] = LeeNMF(MySTFTabs, opts);
apd_muross_time=toc;
apd_muross = W*H;

fprintf('      APD-MURoss Done in time: %f \n',apd_muross_time);

%% EVALUATE APD's MU (Ross)
%inverse transform and cut to size
apd_muross_newsig = istft_catbox(apd_muross.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_muross_newsig = apd_muross_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_muross_W_neg = min(min(W))<0;
apd_muross_H_neg = min(min(H))<0;
fprintf('      APD-MURoss spectra w/ neg values?: %d \n',apd_muross_W_neg);
fprintf('      APD-MURoss coeffs w/ neg values?: %d \n',apd_muross_H_neg);

%compute reconstruction error
apd_muross_rec_err = norm(origMix-apd_muross_newsig)/norm(origMix);
fprintf('      APD-MURoss Normalized Reconstruction Error: %e \n',apd_muross_rec_err);

%compute log spect distance
apd_muross_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_muross(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-MURoss Final Log Spectral Distance: %e \n',apd_muross_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_muross_SDR, apd_muross_SIR, apd_muross_SAR] = bss_eval_sources(apd_muross_newsig',origMix');
fprintf('      APD-MURoss SDR: %f \t SIR: %f \t SAR: %f \n',apd_muross_SDR, apd_muross_SIR, apd_muross_SAR);

disp('--- Finished ---')


%% TEST APD's MU (Brunet) 
% Note: this code seems unadapted
% function [w,h]=nmf(v,r,verbose)

tic;
[ W, H ] = MUNMF(MySTFTabs, nrcomponents, myverbose);
apd_mubrunet_time=toc;
apd_mubrunet = W*H;

fprintf('      APD-MUBrunet Done in time: %f \n',apd_mubrunet_time);

%% EVALUATE APD's MU (Brunet)
%inverse transform and cut to size
apd_mubrunet_newsig = istft_catbox(apd_mubrunet.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_mubrunet_newsig = apd_mubrunet_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_mubrunet_W_neg = min(min(W))<0;
apd_mubrunet_H_neg = min(min(H))<0;
fprintf('      APD-MUBrunet spectra w/ neg values?: %d \n',apd_mubrunet_W_neg);
fprintf('      APD-MUBrunet coeffs w/ neg values?: %d \n',apd_mubrunet_H_neg);

%compute reconstruction error
apd_mubrunet_rec_err = norm(origMix-apd_mubrunet_newsig)/norm(origMix);
fprintf('      APD-MUBrunet Normalized Reconstruction Error: %e \n',apd_mubrunet_rec_err);

%compute log spect distance
apd_mubrunet_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_mubrunet(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-MUBrunet Final Log Spectral Distance: %e \n',apd_mubrunet_logspectdist);


% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_mubrunet_SDR, apd_mubrunet_SIR, apd_mubrunet_SAR] = bss_eval_sources(apd_mubrunet_newsig',origMix');
fprintf('      APD-MUBrunet SDR: %f \t SIR: %f \t SAR: %f \n',apd_mubrunet_SDR, apd_mubrunet_SIR, apd_mubrunet_SAR);

disp('--- Finished ---')

%% NeNMF: also included in lrslibrary with different function header
% skip and use lrslibrary's version as default

%% TEST APD's NtNMF
%function [B, C, HIS] = NtNMF(A, opts, varargin)
% https://www.cs.utexas.edu/~dmkim/Source/software/nnma/nnma.html
% A is the target matrix 
% k is rank of B and C.
% maxit is the maximum number of iterations
% B0 is an initial matrix.
% option = {0|1|2}, 0 for Newton, 1 for Quasi-Newton, 
% 2 for CG. Default is 1.

% Note: tolerance seems to be built in at 1e-2!

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.verbose = myverbose; 
opts.maxIter = maxiter/5;       % note: this is a slow but fast converging impl.

tic;
[ W, H, HIS ] = NtNMF(MySTFTabs, opts);
apd_ntnmf_time=toc;
apd_ntnmf = W*H;

fprintf('      APD-NtNMF Done in time: %f \n',apd_ntnmf_time);

%% EVALUATE APD's NtNMF
%inverse transform and cut to size
apd_ntnmf_newsig = istft_catbox(apd_ntnmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_ntnmf_newsig = apd_ntnmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_ntnmf_W_neg = min(min(W))<0;
apd_ntnmf_H_neg = min(min(H))<0;
fprintf('      APD-NtNMF spectra w/ neg values?: %d \n',apd_ntnmf_W_neg);
fprintf('      APD-NtNMF coeffs w/ neg values?: %d \n',apd_ntnmf_H_neg);

%compute reconstruction error
apd_ntnmf_rec_err = norm(origMix-apd_ntnmf_newsig)/norm(origMix);
fprintf('      APD-NtNMF Normalized Reconstruction Error: %e \n',apd_ntnmf_rec_err);

%compute log spect distance
apd_ntnmf_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_ntnmf(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-NtNMF Final Log Spectral Distance: %e \n',apd_ntnmf_logspectdist);


% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_ntnmf_SDR, apd_ntnmf_SIR, apd_ntnmf_SAR] = bss_eval_sources(apd_ntnmf_newsig',origMix');
fprintf('      APD-NtNMF SDR: %f \t SIR: %f \t SAR: %f \n',apd_ntnmf_SDR, apd_ntnmf_SIR, apd_ntnmf_SAR);

disp('--- Finished ---')


%% TEST APD's PGNMF
% NMF by alternative non-negative least squares using projected gradients
% Source: Lin, Projected Gradient Methods for Nonnegative Matrix Factorization,
% Neural Computation, 19, p. 2756-2779, 2007, MIT press.

% function [W, H, HIS] = PGLIN(V, opts)
% Winit = opts.W0;
% Hinit = opts.H0;
% tol = opts.tolerance;
% timelimit = opts.timeLimit;
% maxiter = opts.maxIter;

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.tolerance = eps;
opts.timeLimit = 10000000;      % set so high that it barely counts
opts.verbose = myverbose; 
opts.maxIter = maxiter;       % note: this is a slow impl.

tic;
[ W, H, HIS ] = PGLIN(MySTFTabs, opts);
apd_pglin_time=toc;
apd_pglin = W*H;

fprintf('      APD-PGLIN Done in time: %f \n',apd_pglin_time);


%% EVALUATE APD's PGNMF
%inverse transform and cut to size
apd_pglin_newsig = istft_catbox(apd_pglin.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_pglin_newsig = apd_pglin_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_pglin_W_neg = min(min(W))<0;
apd_pglin_H_neg = min(min(H))<0;
fprintf('      APD-pglin spectra w/ neg values?: %d \n',apd_pglin_W_neg);
fprintf('      APD-pglin coeffs w/ neg values?: %d \n',apd_pglin_H_neg);

%compute reconstruction error
apd_pglin_rec_err = norm(origMix-apd_pglin_newsig)/norm(origMix);
fprintf('      APD-pglin Normalized Reconstruction Error: %e \n',apd_pglin_rec_err);

%compute log spect distance
apd_pglin_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_pglin(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-pglin Final Log Spectral Distance: %e \n',apd_pglin_logspectdist);


% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_pglin_SDR, apd_pglin_SIR, apd_pglin_SAR] = bss_eval_sources(apd_pglin_newsig',origMix');
fprintf('      APD-pglin SDR: %f \t SIR: %f \t SAR: %f \n',apd_pglin_SDR, apd_pglin_SIR, apd_pglin_SAR);

disp('--- Finished ---')


%% TEST APD's QnNMF
%function [B, C, HIS] = QnNMF(A, opts, varargin)
% https://www.cs.utexas.edu/~dmkim/Source/software/nnma/nnma.html
% A is the target matrix 
% k is rank of B and C.
% maxit is the maximum number of iterations
% B0 is an initial matrix.
% option = {0|1|2}, 0 for Newton, 1 for Quasi-Newton, 
% 2 for CG. Default is 1.

% Note: tolerance seems to be built in at 1e-2!

[m,n] = size(MySTFTabs);
opts.W0 = rand(m,nrcomponents); opts.H0 = rand(nrcomponents,n);
opts.verbose = myverbose; 
opts.maxIter = maxiter/5;       % note: this is a slow impl.      

tic;
[ W, H, HIS ] = QnNMF(MySTFTabs, opts);
apd_qnnmf_time=toc;
apd_qnnmf = W*H;

fprintf('      APD-QnNMF Done in time: %f \n',apd_qnnmf_time);

%% EVALUATE APD's QnNMF
%inverse transform and cut to size
apd_qnnmf_newsig = istft_catbox(apd_qnnmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_qnnmf_newsig = apd_qnnmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_qnnmf_W_neg = min(min(W))<0;
apd_qnnmf_H_neg = min(min(H))<0;
fprintf('      APD-qnnmf spectra w/ neg values?: %d \n',apd_qnnmf_W_neg);
fprintf('      APD-qnnmf coeffs w/ neg values?: %d \n',apd_qnnmf_H_neg);

%compute reconstruction error
apd_qnnmf_rec_err = norm(origMix-apd_qnnmf_newsig)/norm(origMix);
fprintf('      APD-qnnmf Normalized Reconstruction Error: %e \n',apd_qnnmf_rec_err);

%compute log spect distance
apd_qnnmf_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(apd_qnnmf(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      APD-qnnmf Final Log Spectral Distance: %e \n',apd_qnnmf_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_qnnmf_SDR, apd_qnnmf_SIR, apd_qnnmf_SAR] = bss_eval_sources(apd_qnnmf_newsig',origMix');
fprintf('      APD-qnnmf SDR: %f \t SIR: %f \t SAR: %f \n',apd_qnnmf_SDR, apd_qnnmf_SIR, apd_qnnmf_SAR);

disp('--- Finished ---')

%% present results in a table for easier incorporation in other docs

T = table(...
[apd_muacc_time apd_pglinacc_time apd_halsacc_time apd_actset_time apd_blockp_time apd_aloexact_time apd_fcdm_time apd_muross_time apd_mubrunet_time apd_ntnmf_time apd_pglin_time apd_qnnmf_time ]',...
[apd_muacc_W_neg apd_pglinacc_W_neg apd_halsacc_W_neg apd_actset_W_neg apd_blockp_W_neg apd_aloexact_W_neg apd_fcdm_W_neg apd_muross_W_neg apd_mubrunet_W_neg apd_ntnmf_W_neg apd_pglin_W_neg apd_qnnmf_W_neg]',...
[apd_muacc_H_neg apd_pglinacc_H_neg apd_halsacc_H_neg apd_actset_H_neg apd_blockp_H_neg apd_aloexact_H_neg apd_fcdm_H_neg apd_muross_H_neg apd_mubrunet_H_neg apd_ntnmf_H_neg apd_pglin_H_neg apd_qnnmf_H_neg]',...
[apd_muacc_rec_err apd_pglinacc_rec_err apd_halsacc_rec_err apd_actset_rec_err apd_blockp_rec_err apd_aloexact_rec_err apd_fcdm_rec_err apd_muross_rec_err apd_mubrunet_rec_err apd_ntnmf_rec_err apd_pglin_rec_err apd_qnnmf_rec_err]',...
[apd_muacc_logspectdist apd_pglinacc_logspectdist apd_halsacc_logspectdist apd_actset_logspectdist apd_blockp_logspectdist apd_aloexact_logspectdist apd_fcdm_logspectdist apd_muross_logspectdist apd_mubrunet_logspectdist apd_ntnmf_logspectdist apd_pglin_logspectdist apd_qnnmf_logspectdist]',...
[apd_muacc_SDR apd_pglinacc_SDR apd_halsacc_SDR apd_actset_SDR apd_blockp_SDR apd_aloexact_SDR apd_fcdm_SDR apd_muross_SDR apd_mubrunet_SDR apd_ntnmf_SDR apd_pglin_SDR apd_qnnmf_SDR ]',...
'VariableNames', {'Time', 'NegW', 'NegH', 'RecErr', 'LogSpectDist', 'SDR'}, ...
'RowNames', {'APD Accel MU','APD Accel PGLIN', 'APD Accel HALS', 'APD ActiveSet', 'APD BlockPivot', 'APD AloExact', 'APD FCDM', 'APD MURoss', 'APD MUBrunet', 'APD NewtonNMF', 'APD PGLIN', 'APD QNNMF'}...
);

%% delete unnecessary rows/columns
% like those with negative values
toDelete = T.NegW > 0;
T(toDelete, :) = [];
% now delete NegW and NegH columns
T(:,{'NegW', 'NegH'}) = [];
% display table
T