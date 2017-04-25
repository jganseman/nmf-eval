% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% march 2017

% This demo tests different implementations of NMF, as found in the LRS
% library, with audio data.

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

%% TEST lrslibrary's Deep-Semi-Nmf
% function [ Z, H, dnorm ] = deep_seminmf ( X, layers, varargin )
% Optional arguments : 'z0' 'h0' 'bUpdateH' 'bUpdateLastH' 'maxiter' 'TolFun', 'verbose', 'bUpdateZ', 'cache' ...
tic;
[ W, H, dnorm ] = deep_seminmf ( MySTFTabs, nrcomponents, 'maxiter', maxiter, 'verbose', myverbose, 'TolFun', eps);
%note: default TolFun is 10e-5, what is the effect on logarithmic audio data?
lrs_dsnmf_time=toc;

W = cell2mat(W); 
H = cell2mat(H);
lrs_dsnmf=W*H;
%L = W(:,4)*H(4,:);
%imagesc(L); colormap(jet);

fprintf('      LRS-DSNMF Done in time: %f \n',lrs_dsnmf_time);

%% EVALUATE lrslibrary's Deep-Semi-Nmf
%inverse transform and cut to size
lrs_dsnmf_newsig = istft_catbox(lrs_dsnmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
lrs_dsnmf_newsig = lrs_dsnmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
lrs_dsnmf_W_neg = min(min(W))<0;
lrs_dsnmf_H_neg = min(min(H))<0;
fprintf('      LRS-DSNMF spectra w/ neg values?: %d \n',lrs_dsnmf_W_neg);
fprintf('      LRS-DSNMF coeffs w/ neg values?: %d \n',lrs_dsnmf_H_neg);

%compute reconstruction error
rec_err = norm(origMix-lrs_dsnmf_newsig)/norm(origMix);
fprintf('      LRS-DSNMF Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[lrs_dsnmf_SDR lrs_dsnmf_SIR lrs_dsnmf_SAR] = bss_eval_sources(lrs_dsnmf_newsig',origMix');
fprintf('      LRS-DSNMF SDR: %f \t SIR: %f \t SAR: %f \n',lrs_dsnmf_SDR, lrs_dsnmf_SIR, lrs_dsnmf_SAR);

disp('--- Finished ---')

%% TEST lrslibrary's DRMF: Direct Robust Matrix Factorization (Xiong et al. 2011)
%function [ L, S ] = DRMF(M, K, E, options)
% options: 
%     init: the initial L. default the raw SVD.
%     type: 'E' for element outliers, 'R' for row outliers, 'C' for column outliers. default 'E'.
%     max_iter: max iteration. default 100.
%     epsilon: relative change of objective to stop the iteration. default 1e-4.
%     verbose: show info or not

%%% initialization
% Note: LRS uses the following two lines, instead of default init.
% indeed gives better results when initialized with inexact_alm_rpca.
lambda = 1/sqrt(max(size(MySTFTabs)));
L_rpca = inexact_alm_rpca(MySTFTabs, lambda, 1e-5, 10);     %matrix, lambda, tol, maxiter
%sv = svdex(L_rpca);
%rk = EffRank(sv, 0.999);
options.max_iter=maxiter;
options.epsilon=eps;
options.verbose=myverbose;

%%% run
options.init = L_rpca;
tic;
[lrs_drmf,S] = DRMF(MySTFTabs, nrcomponents, 0.1, options);  % 0.1 = max 10% of outliers
% S = full(S);
%%% end
lrs_drmf_time=toc;

% L is our reconstructed low-rank approximation. 
% To obtain individual components, do SVD:
[W, s, H] = svdex(lrs_drmf, nrcomponents);
H=H';

fprintf('      LRS-DRMF Done in time: %f \n',lrs_dsnmf_time);

%% EVALUATE lrslibrary's DRMF
%inverse transform and cut to size
lrs_drmf_newsig = istft_catbox(lrs_drmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
lrs_drmf_newsig = lrs_drmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
lrs_drmf_W_neg = min(min(W))<0;
lrs_drmf_H_neg = min(min(H))<0;
fprintf('      LRS-DRMF spectra w/ neg values?: %d \n',lrs_drmf_W_neg);
fprintf('      LRS-DRMF coeffs w/ neg values?: %d \n',lrs_drmf_H_neg);

%compute reconstruction error
rec_err = norm(origMix-lrs_drmf_newsig)/norm(origMix);
fprintf('      LRS-DRMF Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[lrs_drmf_SDR lrs_drmf_SIR lrs_drmf_SAR] = bss_eval_sources(lrs_drmf_newsig',origMix');
fprintf('      LRS-DRMF SDR: %f \t SIR: %f \t SAR: %f \n',lrs_drmf_SDR, lrs_drmf_SIR, lrs_drmf_SAR);

disp('--- Finished ---')


%% TEST lrslibrary's ENMF: Exact NMF (Gillis and Glineur, 2012)
%function [V,W] = ExactNMF(S, r, max_attempts)

tic;
[W,H] = ExactNMF(MySTFTabs, nrcomponents, 20); % nr attempts has nothing to do with max_iter  
lrs_enmf_time=toc;

lrs_enmf=W*H;

fprintf('      LRS-ENMF Done in time: %f \n',lrs_enmf_time);

%% EVALUATE lrslibrary's ENMF
%inverse transform and cut to size
lrs_enmf_newsig = istft_catbox(lrs_enmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
lrs_enmf_newsig = lrs_enmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
lrs_enmf_W_neg = min(min(W))<0;
lrs_enmf_H_neg = min(min(H))<0;
fprintf('      LRS-ENMF spectra w/ neg values?: %d \n',lrs_enmf_W_neg);
fprintf('      LRS-ENMF coeffs w/ neg values?: %d \n',lrs_enmf_H_neg);

%compute reconstruction error
rec_err = norm(origMix-lrs_enmf_newsig)/norm(origMix);
fprintf('      LRS-ENMF Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[lrs_enmf_SDR lrs_enmf_SIR lrs_enmf_SAR] = bss_eval_sources(lrs_enmf_newsig',origMix');
fprintf('      LRS-ENMF SDR: %f \t SIR: %f \t SAR: %f \n',lrs_enmf_SDR, lrs_enmf_SIR, lrs_enmf_SAR);

disp('--- Finished ---')

%% TEST lrslibrary's iNMF: Incremental Subspace Learning via NMF (Bucak and Gunsel, 2009)
%function [W_new, h, A, B] = inmf( V_new, W, h, A, B, rdim, beta, alpha, maxiter)
%INPUTS:
%   V_new: new data sample, a column vector (d x 1)
%   W: Mixing Matrix -or matrix of basis vectors- (d x rdim)
%   h: initialization for the new encoding vector, a column vector (rdim x 1),
%   A: Matrix to store cummulative V*H' 
%   B: Matrix to store cummulative H*H'
%   rdim: factorization rank
%   beta: weighting coefficient for contribution of initial samples (i.e. A=beta*A+alpha*V_new*h';B=beta*B+alpha*h*h')
%   alpha: weighting coefficient for contribution of the new sample
%   maxiter: maximum number of iterations
%OUTPUTS
%   W_new: Updated mixing matrix 
%   h: Updated encoding vector
%   A, B: Updated A,B matrices

% Note: doing this as in the run_alg.m file in its subdir. It relies on 
% NMF by Hoyer's implementation. Copied in this dir under different name
% to resolve namespace clashes

M = MySTFTabs;

% Execute NMF for the first n samples
n = 10; % take first 10 samples
[W,H] = nmf_boyer(M(:,1:n), nrcomponents, 0, maxiter);
L = W*H;
% Now we can execute iNMF on each new samples
%maxiter = 50;
tic;
A = M(:,1:n)*H';
B = H*H';
h = H(:,end); % Warm start for h
for i = n+1:size(M,2)
  %disp(i);
  M_new = M(:,i);
  [W_new,h,A,B] = inmf(M_new,W,h,A,B,nrcomponents,0.9,0.1,maxiter);
  % H_store(:,i-n) = h; % Just for demonstration
  L(:,end+1) = W_new*h;
end
lrs_inmf_time=toc;
S = M - L;

fprintf('      LRS-iNMF Done in time: %f \n',lrs_inmf_time);

%% EVALUATE lrslibrary's iNMF
lrs_inmf=L;

%inverse transform and cut to size
lrs_inmf_newsig = istft_catbox(lrs_inmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
lrs_inmf_newsig = lrs_inmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
lrs_inmf_W_neg = min(min(W_new))<0;
lrs_inmf_H_neg = min(min(h))<0;
fprintf('      LRS-iNMF spectra w/ neg values?: %d \n',lrs_inmf_W_neg);
fprintf('      LRS-iNMF coeffs w/ neg values?: %d \n',lrs_inmf_H_neg);

%compute reconstruction error
rec_err = norm(origMix-lrs_inmf_newsig)/norm(origMix);
fprintf('      LRS-iNMF Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[lrs_inmf_SDR lrs_inmf_SIR lrs_inmf_SAR] = bss_eval_sources(lrs_inmf_newsig',origMix');
fprintf('      LRS-iNMF SDR: %f \t SIR: %f \t SAR: %f \n',lrs_inmf_SDR, lrs_inmf_SIR, lrs_inmf_SAR);

disp('--- Finished ---')

%this routine made a few extra directories. Force remove them
rmdir('initials', 's');

%% TEST lrslibrary's LMNF: Spatially Localized NMF (Li et al. 2001)
% Reference:
%   S. Z. Li, X. Hou, H. Zhang, and Q. Cheng, "Learning spatially
%   localized, parts-based representation," in IEEE Conference on Computer Vision and Pattern Recognition, 2001, pp. 207–212.

% Relies on a Mex-function to compute kullback-leibler divergence: KLC.cpp
% install a C++ compiler and run mex -setup 
% then compile (run mex) in /lrslibrary/algorithms/nmf/LNMF
% (function signature was erroneous so fixed and copied in this dir)

option.verbose = myverbose;
% implementation from LRS has the following options built-in:
% niter = 3000;     % maximum number of iterations
% precision = 1e-4;
tic;
[W,H] = LNMF(MySTFTabs,nrcomponents,option);
lrs_lnmf_time=toc;

lrs_lnmf = W * H;
S = MySTFTabs - lrs_lnmf;

fprintf('      LRS-LNMF Done in time: %f \n',lrs_lnmf_time);

%% EVALUATE lrslibrary's LNMF

%inverse transform and cut to size
lrs_lnmf_newsig = istft_catbox(lrs_lnmf.*phase, fftsize / hopsize, fftsize, 'smooth')';
lrs_lnmf_newsig = lrs_lnmf_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
lrs_lnmf_W_neg = min(min(W))<0;
lrs_lnmf_H_neg = min(min(H))<0;
fprintf('      LRS-LNMF spectra w/ neg values?: %d \n',lrs_lnmf_W_neg);
fprintf('      LRS-LNMF coeffs w/ neg values?: %d \n',lrs_lnmf_H_neg);

%compute reconstruction error
rec_err = norm(origMix-lrs_lnmf_newsig)/norm(origMix);
fprintf('      LRS-LNMF Normalized Reconstruction Error: %e \n',rec_err);

% compute BSS EVAL. make sure we define rows as signals, not columns
[lrs_lnmf_SDR lrs_lnmf_SIR lrs_lnmf_SAR] = bss_eval_sources(lrs_lnmf_newsig',origMix');
fprintf('      LRS-iNMF SDR: %f \t SIR: %f \t SAR: %f \n',lrs_lnmf_SDR, lrs_lnmf_SIR, lrs_lnmf_SAR);

disp('--- Finished ---')
