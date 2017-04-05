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