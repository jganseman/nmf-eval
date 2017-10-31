% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% october 2017

% This demo tests different implementations of NMF, as found in Graham
% Grindlay's nmflib , with audio data.

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
maxiter = 50;
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

%% feed algos with same initial matrices
[numrows,numcols] = size(MySTFTabs);
MyW0 = rand(numrows,nrcomponents);
MyH0 = rand(nrcomponents, numcols);

%% TEST nmflib's nmf_amari
% function [W,H,errs,vout] = nmf_amari(V,r,varargin)
tic;
[ W, H, gri_amari_errs ] = nmf_amari ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'alpha', 0.5);
gri_amari_time=toc;
gri_amari=W*H;
fprintf('      GRI-AMARI @ 0.5 Done in time: %f \n',gri_amari_time);

%% EVALUATE nmflib's nmf_amari
%inverse transform and cut to size
gri_amari_newsig = istft_catbox(gri_amari.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_amari_newsig = gri_amari_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_amari_W_neg = min(min(W))<0;
gri_amari_H_neg = min(min(H))<0;
fprintf('      GRI-AMARI @ 0.5 spectra w/ neg values?: %d \n',gri_amari_W_neg);
fprintf('      GRI-AMARI @ 0.5 coeffs w/ neg values?: %d \n',gri_amari_H_neg);

%compute reconstruction error
gri_amari_rec_err = norm(origMix-gri_amari_newsig)/norm(origMix);
fprintf('      GRI-AMARI @ 0.5 Normalized Reconstruction Error: %e \n',gri_amari_rec_err);

%compute log spect distance
gri_amari_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_amari(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-AMARI @ 0.5 Final Log Spectral Distance: %e \n',gri_amari_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_amari_SDR, gri_amari_SIR, gri_amari_SAR] = bss_eval_sources(gri_amari_newsig',origMix');
fprintf('      GRI-AMARI @ 0.5 SDR: %f \t SIR: %f \t SAR: %f \n',gri_amari_SDR, gri_amari_SIR, gri_amari_SAR);
disp('--- Finished ---')


%% TEST nmflib's nmf_beta
% function [W,H,errs,vout] = nmf_beta(V,r,varargin)
tic;
[ W, H, gri_beta_errs ] = nmf_beta ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'beta', 1.0);
gri_beta_time=toc;
gri_beta=W*H;
fprintf('      GRI-BETA @ 1.0 Done in time: %f \n',gri_beta_time);

%% EVALUATE nmflib's nmf_beta
%inverse transform and cut to size
gri_beta_newsig = istft_catbox(gri_beta.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_beta_newsig = gri_beta_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_beta_W_neg = min(min(W))<0;
gri_beta_H_neg = min(min(H))<0;
fprintf('      GRI-BETA @ 1.0 spectra w/ neg values?: %d \n',gri_beta_W_neg);
fprintf('      GRI-BETA @ 1.0 coeffs w/ neg values?: %d \n',gri_beta_H_neg);

%compute reconstruction error
gri_beta_rec_err = norm(origMix-gri_beta_newsig)/norm(origMix);
fprintf('      GRI-BETA @ 1.0 Normalized Reconstruction Error: %e \n',gri_beta_rec_err);

%compute log spect distance
gri_beta_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_beta(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-BETA @ 1.0 Final Log Spectral Distance: %e \n',gri_beta_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_beta_SDR, gri_beta_SIR, gri_beta_SAR] = bss_eval_sources(gri_beta_newsig',origMix');
fprintf('      GRI-BETA @ 1.0 SDR: %f \t SIR: %f \t SAR: %f \n',gri_beta_SDR, gri_beta_SIR, gri_beta_SAR);
disp('--- Finished ---')

%% TEST nmflib's nmf_convex
% function [W,H,errs,vout] = nmf_convex(V, r, varargin)
% tic;
% [ W, H, gri_convex_errs ] = nmf_convex ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0);
% gri_convex_time=toc;
% gri_convex=W*H;
% fprintf('      GRI-CONVEX Done in time: %f \n',gri_convex_time);

%NOTE: does not work! Errors in nmflib code...
%But can give negative results anyway, so don't do this method?

%% EVALUATE nmflib's nmf_convex

%% TEST nmflib's nmf_euc
tic;
[ W, H, gri_euc_errs ] = nmf_euc ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0);
gri_euc_time=toc;
gri_euc=W*H;
fprintf('      GRI-EUC Done in time: %f \n',gri_euc_time);

%% EVALUATE nmflib's nmf_euc
%inverse transform and cut to size
gri_euc_newsig = istft_catbox(gri_euc.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_euc_newsig = gri_euc_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_euc_W_neg = min(min(W))<0;
gri_euc_H_neg = min(min(H))<0;
fprintf('      GRI-EUC spectra w/ neg values?: %d \n',gri_euc_W_neg);
fprintf('      GRI-EUC coeffs w/ neg values?: %d \n',gri_euc_H_neg);

%compute reconstruction error
gri_euc_rec_err = norm(origMix-gri_euc_newsig)/norm(origMix);
fprintf('      GRI-EUC Normalized Reconstruction Error: %e \n',gri_euc_rec_err);

%compute log spect distance
gri_euc_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_euc(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-EUC Final Log Spectral Distance: %e \n',gri_euc_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_euc_SDR, gri_euc_SIR, gri_euc_SAR] = bss_eval_sources(gri_euc_newsig',origMix');
fprintf('      GRI-EUC SDR: %f \t SIR: %f \t SAR: %f \n',gri_euc_SDR, gri_euc_SIR, gri_euc_SAR);
disp('--- Finished ---')


%% TEST nmflib's nmf_euc_orth
% accepts bools orth_w (default true) and orth_h (default false)
myorthw = 1;
myorthh = 0;
tic;
[ W, H, gri_orth_errs ] = nmf_euc_orth ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'orth_w', myorthw, 'orth_h', myorthh);
gri_orth_time=toc;
gri_orth=W*H;
fprintf('      GRI-ORTH Done in time: %f \n',gri_orth_time);

%% EVALUATE nmflib's nmf_euc_orth
%inverse transform and cut to size
gri_orth_newsig = istft_catbox(gri_orth.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_orth_newsig = gri_orth_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_orth_W_neg = min(min(W))<0;
gri_orth_H_neg = min(min(H))<0;
fprintf('      GRI-ORTH spectra w/ neg values?: %d \n',gri_orth_W_neg);
fprintf('      GRI-ORTH coeffs w/ neg values?: %d \n',gri_orth_H_neg);

%compute reconstruction error
gri_orth_rec_err = norm(origMix-gri_orth_newsig)/norm(origMix);
fprintf('      GRI-ORTH Normalized Reconstruction Error: %e \n',gri_orth_rec_err);

%compute log spect distance
gri_orth_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_orth(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-ORTH Final Log Spectral Distance: %e \n',gri_orth_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_orth_SDR, gri_orth_SIR, gri_orth_SAR] = bss_eval_sources(gri_orth_newsig',origMix');
fprintf('      GRI-ORTH SDR: %f \t SIR: %f \t SAR: %f \n',gri_orth_SDR, gri_orth_SIR, gri_orth_SAR);
disp('--- Finished ---')


%% TEST nmflib's nmf_euc_sparse_es
% function [W,H,errs,vout] = nmf_euc_sparse_es(V, r, varargin)
% default sparsity parameter: alpha = 0
myalpha = 0.1;
tic;
[ W, H, gri_eucsp_errs ] = nmf_euc_sparse_es ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'alpha', myalpha);
gri_eucsp_time=toc;
gri_eucsp=W*H;
fprintf('      GRI-EUCSP Done in time: %f \n',gri_eucsp_time);

%% EVALUATE nmflib's nmf_euc_sparse_es
%inverse transform and cut to size
gri_eucsp_newsig = istft_catbox(gri_eucsp.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_eucsp_newsig = gri_eucsp_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_eucsp_W_neg = min(min(W))<0;
gri_eucsp_H_neg = min(min(H))<0;
fprintf('      GRI-EUCSP spectra w/ neg values?: %d \n',gri_eucsp_W_neg);
fprintf('      GRI-EUCSP coeffs w/ neg values?: %d \n',gri_eucsp_H_neg);

%compute reconstruction error
gri_eucsp_rec_err = norm(origMix-gri_eucsp_newsig)/norm(origMix);
fprintf('      GRI-EUCSP Normalized Reconstruction Error: %e \n',gri_eucsp_rec_err);

%compute log spect distance
gri_eucsp_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_eucsp(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-EUCSP Final Log Spectral Distance: %e \n',gri_eucsp_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_eucsp_SDR, gri_eucsp_SIR, gri_eucsp_SAR] = bss_eval_sources(gri_eucsp_newsig',origMix');
fprintf('      GRI-EUCSP SDR: %f \t SIR: %f \t SAR: %f \n',gri_eucsp_SDR, gri_eucsp_SIR, gri_eucsp_SAR);
disp('--- Finished ---')


%% TEST nmflib's nmf_kl
tic;
[ W, H, gri_kl_errs ] = nmf_kl ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0);
gri_kl_time=toc;
gri_kl=W*H;
fprintf('      GRI-KL Done in time: %f \n',gri_kl_time);

%% EVALUATE nmflib's nmf_kl
%inverse transform and cut to size
gri_kl_newsig = istft_catbox(gri_kl.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_kl_newsig = gri_kl_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_kl_W_neg = min(min(W))<0;
gri_kl_H_neg = min(min(H))<0;
fprintf('      GRI-KL spectra w/ neg values?: %d \n',gri_kl_W_neg);
fprintf('      GRI-KL coeffs w/ neg values?: %d \n',gri_kl_H_neg);

%compute reconstruction error
gri_kl_rec_err = norm(origMix-gri_kl_newsig)/norm(origMix);
fprintf('      GRI-KL Normalized Reconstruction Error: %e \n',gri_kl_rec_err);

%compute log spect distance
gri_kl_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_kl(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-KL Final Log Spectral Distance: %e \n',gri_kl_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_kl_SDR, gri_kl_SIR, gri_kl_SAR] = bss_eval_sources(gri_kl_newsig',origMix');
fprintf('      GRI-KL SDR: %f \t SIR: %f \t SAR: %f \n',gri_kl_SDR, gri_kl_SIR, gri_kl_SAR);
disp('--- Finished ---')

%% TEST nmflib's nmf_kl_con
% convolutive NMF
% extra parameter: 'win' = width of columns W, default 1.
% returns matrix W as (n x r x win)
% to reconstruct result, see function R = rec_cnmf(W,H,myeps)
mywin = 1;
tic;
[ W, H, gri_klcon_errs ] = nmf_kl_con ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'win', mywin);
gri_klcon_time=toc;
gri_klcon=W*H;
fprintf('      GRI-KL-CON Done in time: %f \n',gri_klcon_time);
% Note: currently does not work with windows larger than 10

%% EVALUATE nmflib's nmf_kl_con
gri_klcon_newsig = istft_catbox(gri_klcon.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_klcon_newsig = gri_klcon_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_klcon_W_neg = min(min(W))<0;
gri_klcon_H_neg = min(min(H))<0;
fprintf('      GRI-KL-CON spectra w/ neg values?: %d \n',gri_klcon_W_neg);
fprintf('      GRI-KL-CON coeffs w/ neg values?: %d \n',gri_klcon_H_neg);

%compute reconstruction error
gri_klcon_rec_err = norm(origMix-gri_klcon_newsig)/norm(origMix);
fprintf('      GRI-KL-CON Normalized Reconstruction Error: %e \n',gri_klcon_rec_err);

%compute log spect distance
gri_klcon_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_kl(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-KL-CON Final Log Spectral Distance: %e \n',gri_klcon_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_klcon_SDR, gri_klcon_SIR, gri_klcon_SAR] = bss_eval_sources(gri_klcon_newsig',origMix');
fprintf('      GRI-KL-CON SDR: %f \t SIR: %f \t SAR: %f \n',gri_klcon_SDR, gri_klcon_SIR, gri_klcon_SAR);
disp('--- Finished ---')

%% TEST nmflib's nmf_kl_loc
% sparse (spatially localized) version with constraints according to S. Li (2001)
tic;
[ W, H, gri_klloc_errs ] = nmf_kl_loc ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0);
gri_klloc_time=toc;
gri_klloc=W*H;
fprintf('      GRI-KL-LOC Done in time: %f \n',gri_klloc_time);

%% EVALUATE nmflib's nmf_kl_loc
%inverse transform and cut to size
gri_klloc_newsig = istft_catbox(gri_klloc.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_klloc_newsig = gri_klloc_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_klloc_W_neg = min(min(W))<0;
gri_klloc_H_neg = min(min(H))<0;
fprintf('      GRI-KL-LOC spectra w/ neg values?: %d \n',gri_klloc_W_neg);
fprintf('      GRI-KL-LOC coeffs w/ neg values?: %d \n',gri_klloc_H_neg);

%compute reconstruction error
gri_klloc_rec_err = norm(origMix-gri_klloc_newsig)/norm(origMix);
fprintf('      GRI-KL-LOC Normalized Reconstruction Error: %e \n',gri_klloc_rec_err);

%compute log spect distance
gri_klloc_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_klloc(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-KL-LOC Final Log Spectral Distance: %e \n',gri_klloc_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_klloc_SDR, gri_klloc_SIR, gri_klloc_SAR] = bss_eval_sources(gri_klloc_newsig',origMix');
fprintf('      GRI-KL-LOC SDR: %f \t SIR: %f \t SAR: %f \n',gri_klloc_SDR, gri_klloc_SIR, gri_klloc_SAR);
disp('--- Finished ---')


%% TEST nmflib's nmf_kl_ns
% does a tri-factorization, V = WSH, with S defined by a smoothing parameter alpha, default 0
mysmooth = 0.1;
tic;
[ W, H, gri_klns_errs ] = nmf_kl_ns ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'alpha', mysmooth);
gri_klns_time=toc;
gri_klns=W*H;
fprintf('      GRI-KL-NS Done in time: %f \n',gri_klns_time);

%% EVALUATE nmflib's nmf_kl_ns
%inverse transform and cut to size
gri_klns_newsig = istft_catbox(gri_klns.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_klns_newsig = gri_klns_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_klns_W_neg = min(min(W))<0;
gri_klns_H_neg = min(min(H))<0;
fprintf('      GRI-KL-NS spectra w/ neg values?: %d \n',gri_klns_W_neg);
fprintf('      GRI-KL-NS coeffs w/ neg values?: %d \n',gri_klns_H_neg);

%compute reconstruction error
gri_klns_rec_err = norm(origMix-gri_klns_newsig)/norm(origMix);
fprintf('      GRI-KL-NS Normalized Reconstruction Error: %e \n',gri_klns_rec_err);

%compute log spect distance
gri_klns_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_klns(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-KL-NS Final Log Spectral Distance: %e \n',gri_klns_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_klns_SDR, gri_klns_SIR, gri_klns_SAR] = bss_eval_sources(gri_klns_newsig',origMix');
fprintf('      GRI-KL-NS SDR: %f \t SIR: %f \t SAR: %f \n',gri_klns_SDR, gri_klns_SIR, gri_klns_SAR);
disp('--- Finished ---')

%% TEST nmflib's nmf_kl_sparse_es
% sparsity according to Eggert / Schmidt constraints
% uses parameter alpha to control sparsity, defaults to 0
tic;
[ W, H, gri_klsp_errs ] = nmf_kl_sparse_es ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'alpha', myalpha);
gri_klsp_time=toc;
gri_klsp=W*H;
fprintf('      GRI-KL-SP Done in time: %f \n',gri_klsp_time);

%% EVALUATE nmflib's nmf_kl_sparse_es
%inverse transform and cut to size
gri_klsp_newsig = istft_catbox(gri_klsp.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_klsp_newsig = gri_klsp_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_klsp_W_neg = min(min(W))<0;
gri_klsp_H_neg = min(min(H))<0;
fprintf('      GRI-KL-SP spectra w/ neg values?: %d \n',gri_klsp_W_neg);
fprintf('      GRI-KL-SP coeffs w/ neg values?: %d \n',gri_klsp_H_neg);

%compute reconstruction error
gri_klsp_rec_err = norm(origMix-gri_klsp_newsig)/norm(origMix);
fprintf('      GRI-KL-SP Normalized Reconstruction Error: %e \n',gri_klsp_rec_err);

%compute log spect distance
gri_klsp_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_klsp(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-KL-SP Final Log Spectral Distance: %e \n',gri_klsp_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_klsp_SDR, gri_klsp_SIR, gri_klsp_SAR] = bss_eval_sources(gri_klsp_newsig',origMix');
fprintf('      GRI-KL-SP SDR: %f \t SIR: %f \t SAR: %f \n',gri_klsp_SDR, gri_klsp_SIR, gri_klsp_SAR);
disp('--- Finished ---')

%% TEST nmflib's nmf_kl_sparse_v
% sparsity according to Virtanen constraints
% uses parameter alpha to control sparsity, defaults to 0
tic;
[ W, H, gri_klspv_errs ] = nmf_kl_sparse_v ( MySTFTabs, nrcomponents, 'niter', maxiter, 'verb', myverbose, 'W0', MyW0, 'H0', MyH0, 'alpha', myalpha);
gri_klspv_time=toc;
gri_klspv=W*H;
fprintf('      GRI-KL-SPV Done in time: %f \n',gri_klspv_time);

%% EVALUATE nmflib's nmf_kl_sparse_v

%inverse transform and cut to size
gri_klspv_newsig = istft_catbox(gri_klspv.*phase, fftsize / hopsize, fftsize, 'smooth')';
gri_klspv_newsig = gri_klspv_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
gri_klspv_W_neg = min(min(W))<0;
gri_klspv_H_neg = min(min(H))<0;
fprintf('      GRI-KL-SPV spectra w/ neg values?: %d \n',gri_klspv_W_neg);
fprintf('      GRI-KL-SPV coeffs w/ neg values?: %d \n',gri_klspv_H_neg);

%compute reconstruction error
gri_klspv_rec_err = norm(origMix-gri_klspv_newsig)/norm(origMix);
fprintf('      GRI-KL-SPV Normalized Reconstruction Error: %e \n',gri_klspv_rec_err);

%compute log spect distance
gri_klspv_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(gri_klspv(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      GRI-KL-SPV Final Log Spectral Distance: %e \n',gri_klspv_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[gri_klspv_SDR, gri_klspv_SIR, gri_klspv_SAR] = bss_eval_sources(gri_klspv_newsig',origMix');
fprintf('      GRI-KL-SPV SDR: %f \t SIR: %f \t SAR: %f \n',gri_klspv_SDR, gri_klspv_SIR, gri_klspv_SAR);
disp('--- Finished ---')



%% present results in a table for easier incorporation in other docs

T = table(...
[gri_amari_time gri_beta_time gri_euc_time gri_orth_time gri_eucsp_time gri_kl_time gri_klcon_time gri_klloc_time gri_klns_time gri_klsp_time gri_klspv_time]',...
[gri_amari_W_neg gri_beta_W_neg gri_euc_W_neg gri_orth_W_neg gri_eucsp_W_neg gri_kl_W_neg gri_klcon_W_neg gri_klloc_W_neg gri_klns_W_neg gri_klsp_W_neg gri_klspv_W_neg]',...
[gri_amari_H_neg gri_beta_H_neg gri_euc_H_neg gri_orth_H_neg gri_eucsp_H_neg gri_kl_H_neg gri_klcon_H_neg gri_klloc_H_neg gri_klns_H_neg gri_klsp_H_neg gri_klspv_H_neg]',...
[gri_amari_rec_err gri_beta_rec_err gri_euc_rec_err gri_orth_rec_err gri_eucsp_rec_err gri_kl_rec_err gri_klcon_rec_err gri_klloc_rec_err gri_klns_rec_err gri_klsp_rec_err gri_klspv_rec_err]',...
[gri_amari_logspectdist gri_beta_logspectdist gri_euc_logspectdist gri_orth_logspectdist gri_eucsp_logspectdist gri_kl_logspectdist gri_klcon_logspectdist gri_klloc_logspectdist gri_klns_logspectdist gri_klsp_logspectdist gri_klspv_logspectdist]',...
[gri_amari_SDR gri_beta_SDR gri_euc_SDR gri_orth_SDR gri_eucsp_SDR gri_kl_SDR gri_klcon_SDR gri_klloc_SDR gri_klns_SDR gri_klsp_SDR gri_klspv_SDR]',...
'VariableNames', {'Time', 'NegW', 'NegH', 'RecErr', 'LogSpectDist', 'SDR'}, ...
'RowNames', {'gri Amari','gri Beta', 'gri Euc', 'gri EucOrth', 'gri EucSp', 'gri KL', 'gri KL Convol', 'gri KL Local', 'gri KL Nonsmth', 'gri KLSP_ES', 'gri KLSP_V'}...
);
 
%% delete unnecessary rows/columns
% like those with negative values
toDelete = T.NegW > 0;
T(toDelete, :) = [];
% now delete NegW and NegH columns
T(:,{'NegW', 'NegH'}) = [];
% display table
T


%% plot evolution of errors
figure; 
plot([gri_amari_errs gri_beta_errs] );
legend({'nmflib Amari','nmflib Beta'});
figure;
plot([gri_euc_errs gri_orth_errs gri_eucsp_errs]);
legend({'nmflib Euc', 'nmflib EucOrth', 'nmflib EucSp'});
figure; 
%don't plot convolutive (gri_klcon_errs) with window 1, is same as KL 
plot([gri_kl_errs gri_klloc_errs gri_klns_errs gri_klsp_errs gri_klspv_errs]);
legend({'nmflib KL', 'nmflib KL Local', 'nmflib KL Nonsmth', 'nmflib KL SP-ES', 'nmflib KL SP-V'});

