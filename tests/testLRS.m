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
addpath(genpath('../'));     % add subdirectories to path


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

