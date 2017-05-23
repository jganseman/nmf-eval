% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% may 2017

% This demo tests different implementations of NMF, as found in the 
% NMF_APD library, with audio data.

% For the STFT, we choose the corrected version from CATbox, which is available
% in github.com/jganseman/nsgt_eval . This is the fastest STFT among the 4 
% we tested, and also the one with the smallest reconstruction error when
% window-corrected (option 'smooth').

% Phase is ignored in the NMF; the inverse STFT uses the original phase data.

% For separation instead of signal reconstruction, BSS_EVAL 2.1 is used.

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

fprintf('      APD_MUacc Done in time: %f \n',apd_muacc_time);

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
rec_err = norm(origMix-apd_muacc_newsig)/norm(origMix);
fprintf('      APD-MUacc Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_muacc_SDR apd_muacc_SIR apd_muacc_SAR] = bss_eval_sources(apd_muacc_newsig',origMix');
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

fprintf('      APD_PGLINacc Done in time: %f \n',apd_pglinacc_time);

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
rec_err = norm(origMix-apd_pglinacc_newsig)/norm(origMix);
fprintf('      APD-PGLINacc Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_pglinacc_SDR apd_pglinacc_SIR apd_pglinacc_SAR] = bss_eval_sources(apd_pglinacc_newsig',origMix');
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

fprintf('      APD_HALSacc Done in time: %f \n',apd_halsacc_time);

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
rec_err = norm(origMix-apd_halsacc_newsig)/norm(origMix);
fprintf('      APD-HALSacc Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_halsacc_SDR apd_halsacc_SIR apd_halsacc_SAR] = bss_eval_sources(apd_halsacc_newsig',origMix');
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

fprintf('      APD_ActiveSet Done in time: %f \n',apd_actset_time);

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
rec_err = norm(origMix-apd_actset_newsig)/norm(origMix);
fprintf('      APD-ActiveSet Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_actset_SDR apd_actset_SIR apd_actset_SAR] = bss_eval_sources(apd_actset_newsig',origMix');
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

fprintf('      APD_BlockPivot Done in time: %f \n',apd_blockp_time);

%% EVALUATE APD's BlockPivot
%inverse transform and cut to size
apd_blockp_newsig = istft_catbox(apd_blockp.*phase, fftsize / hopsize, fftsize, 'smooth')';
apd_blockp_newsig = apd_blockp_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
apd_blockp_W_neg = min(min(W))<0;
apd_blockp_H_neg = min(min(H))<0;
fprintf('      APD-ActiveSet spectra w/ neg values?: %d \n',apd_blockp_W_neg);
fprintf('      APD-ActiveSet coeffs w/ neg values?: %d \n',apd_blockp_H_neg);

%compute reconstruction error
rec_err = norm(origMix-apd_blockp_newsig)/norm(origMix);
fprintf('      APD-ActiveSet Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[apd_blockp_SDR apd_blockp_SIR apd_blockp_SAR] = bss_eval_sources(apd_blockp_newsig',origMix');
fprintf('      APD-ActiveSet SDR: %f \t SIR: %f \t SAR: %f \n',apd_blockp_SDR, apd_blockp_SIR, apd_blockp_SAR);

disp('--- Finished ---')