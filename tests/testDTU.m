% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% april 2017

% This demo tests different implementations of NMF, as found in the 
% NMF_DTU library, with audio data.

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

%% TEST DTU's nmf_als
% function [W,H]=nmf_als(X,K,Nitsmax,speak)

tic;
[ W, H ] = nmf_als ( MySTFTabs, nrcomponents, maxiter, myverbose);
dtu_als_time=toc;
dtu_als = W*H;

fprintf('      DTU-ALS Done in time: %f \n',dtu_als_time);

%% EVALUATE DTU's nmf_als
%inverse transform and cut to size
dtu_als_newsig = istft_catbox(dtu_als.*phase, fftsize / hopsize, fftsize, 'smooth')';
dtu_als_newsig = dtu_als_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
dtu_als_W_neg = min(min(W))<0;
dtu_als_H_neg = min(min(H))<0;
fprintf('      DTU-ALS spectra w/ neg values?: %d \n',dtu_als_W_neg);
fprintf('      DTU-ALS coeffs w/ neg values?: %d \n',dtu_als_H_neg);

%compute reconstruction error
dtu_als_rec_err = norm(origMix-dtu_als_newsig)/norm(origMix);
fprintf('      DTU-ALS Normalized Reconstruction Error: %e \n',dtu_als_rec_err);

%compute log spect distance
dtu_als_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(dtu_als(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      DTU-ALS Final Log Spectral Distance: %e \n',dtu_als_logspectdist);


% compute BSS EVAL. make sure we define rows as signals, not columns
[dtu_als_SDR, dtu_als_SIR, dtu_als_SAR] = bss_eval_sources(dtu_als_newsig',origMix');
fprintf('      DTU-ALS SDR: %f \t SIR: %f \t SAR: %f \n',dtu_als_SDR, dtu_als_SIR, dtu_als_SAR);

disp('--- Finished ---')

%% TEST DTU's nmf_alsobs
% function [W,H]=nmf_alsobs(X,K,maxiter,speak)

tic;
[ W, H ] = nmf_alsobs ( MySTFTabs, nrcomponents, maxiter, myverbose);
dtu_alsobs_time=toc;
dtu_alsobs = W*H;

fprintf('      DTU-ALSOBS Done in time: %f \n',dtu_alsobs_time);

%% EVALUATE DTU's nmf_alsobs
%inverse transform and cut to size
dtu_alsobs_newsig = istft_catbox(dtu_alsobs.*phase, fftsize / hopsize, fftsize, 'smooth')';
dtu_alsobs_newsig = dtu_alsobs_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
dtu_alsobs_W_neg = min(min(W))<0;
dtu_alsobs_H_neg = min(min(H))<0;
fprintf('      DTU-ALSOBS spectra w/ neg values?: %d \n',dtu_alsobs_W_neg);
fprintf('      DTU-ALSOBS coeffs w/ neg values?: %d \n',dtu_alsobs_H_neg);

%compute reconstruction error
dtu_alsobs_rec_err = norm(origMix-dtu_alsobs_newsig)/norm(origMix);
fprintf('      DTU-ALSOBS Normalized Reconstruction Error: %e \n',dtu_alsobs_rec_err);

%compute log spect distance
dtu_alsobs_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(dtu_alsobs(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      DTU-ALSOBS Final Log Spectral Distance: %e \n',dtu_alsobs_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[dtu_alsobs_SDR, dtu_alsobs_SIR, dtu_alsobs_SAR] = bss_eval_sources(dtu_alsobs_newsig',origMix');
fprintf('      DTU-ALSOBS SDR: %f \t SIR: %f \t SAR: %f \n',dtu_alsobs_SDR, dtu_alsobs_SIR, dtu_alsobs_SAR);

disp('--- Finished ---')

%% TEST DTU's nmf_cjlin
% function [W,H] = nmf_cjlin(V,Winit,Hinit,tol,timelimit,maxiter)

[D,N] = size(MySTFTabs);
Winit = rand(D,nrcomponents);
Hinit = rand(nrcomponents, N);
timelimit = 10000;

tic;
[ W, H ] = nmf_cjlin (MySTFTabs, Winit, Hinit, eps, timelimit, maxiter);
dtu_cjlin_time=toc;
dtu_cjlin = W*H;

fprintf('      DTU-CJLIN Done in time: %f \n',dtu_cjlin_time);

%% EVALUATE DTU's nmf_cjlin
%inverse transform and cut to size
dtu_cjlin_newsig = istft_catbox(dtu_cjlin.*phase, fftsize / hopsize, fftsize, 'smooth')';
dtu_cjlin_newsig = dtu_cjlin_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
dtu_cjlin_W_neg = min(min(W))<0;
dtu_cjlin_H_neg = min(min(H))<0;
fprintf('      DTU-CJLIN spectra w/ neg values?: %d \n',dtu_cjlin_W_neg);
fprintf('      DTU-CJLIN coeffs w/ neg values?: %d \n',dtu_cjlin_H_neg);

%compute reconstruction error
dtu_cjlin_rec_err = norm(origMix-dtu_cjlin_newsig)/norm(origMix);
fprintf('      DTU-CJLIN Normalized Reconstruction Error: %e \n',dtu_cjlin_rec_err);

%compute log spect distance
dtu_cjlin_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(dtu_cjlin(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      DTU-CJLIN Final Log Spectral Distance: %e \n',dtu_cjlin_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[dtu_cjlin_SDR, dtu_cjlin_SIR, dtu_cjlin_SAR] = bss_eval_sources(dtu_cjlin_newsig',origMix');
fprintf('      DTU-CJLIN SDR: %f \t SIR: %f \t SAR: %f \n',dtu_cjlin_SDR, dtu_cjlin_SIR, dtu_cjlin_SAR);

disp('--- Finished ---')

%% TEST DTU's nmf_mm
% function [W,H]=nmf_mm(X,K,Nitsmax,speak)

tic;
[ W, H ] = nmf_mm ( MySTFTabs, nrcomponents, maxiter, myverbose);
dtu_mm_time=toc;
dtu_mm = W*H;

fprintf('      DTU-MM Done in time: %f \n',dtu_mm_time);

%% EVALUATE DTU's nmf_mm
%inverse transform and cut to size
dtu_mm_newsig = istft_catbox(dtu_mm.*phase, fftsize / hopsize, fftsize, 'smooth')';
dtu_mm_newsig = dtu_mm_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
dtu_mm_W_neg = min(min(W))<0;
dtu_mm_H_neg = min(min(H))<0;
fprintf('      DTU-MM spectra w/ neg values?: %d \n',dtu_mm_W_neg);
fprintf('      DTU-MM coeffs w/ neg values?: %d \n',dtu_mm_H_neg);

%compute reconstruction error
dtu_mm_rec_err = norm(origMix-dtu_mm_newsig)/norm(origMix);
fprintf('      DTU-MM Normalized Reconstruction Error: %e \n',dtu_mm_rec_err);

%compute log spect distance
dtu_mm_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(dtu_mm(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      DTU-MM Final Log Spectral Distance: %e \n',dtu_mm_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[dtu_mm_SDR, dtu_mm_SIR, dtu_mm_SAR] = bss_eval_sources(dtu_mm_newsig',origMix');
fprintf('      DTU-MM SDR: %f \t SIR: %f \t SAR: %f \n',dtu_mm_SDR, dtu_mm_SIR, dtu_mm_SAR);

disp('--- Finished ---')

%% TEST DTU's nmf_prob
% function [W,H]=nmf_prob(X,K,Nitsmax,speak)

tic;
[ W, H ] = nmf_prob ( MySTFTabs, nrcomponents, maxiter, myverbose);
dtu_prob_time=toc;
dtu_prob = W*H;

fprintf('      DTU-PROB Done in time: %f \n',dtu_prob_time);

%% EVALUATE DTU's nmf_prob
%inverse transform and cut to size
dtu_prob_newsig = istft_catbox(dtu_prob.*phase, fftsize / hopsize, fftsize, 'smooth')';
dtu_prob_newsig = dtu_prob_newsig(fftsize+1:fftsize+length(origMix));

%has negative values?
dtu_prob_W_neg = min(min(W))<0;
dtu_prob_H_neg = min(min(H))<0;
fprintf('      DTU-PROB spectra w/ neg values?: %d \n',dtu_prob_W_neg);
fprintf('      DTU-PROB coeffs w/ neg values?: %d \n',dtu_prob_H_neg);

%compute reconstruction error
dtu_prob_rec_err = norm(origMix-dtu_prob_newsig)/norm(origMix);
fprintf('      DTU-PROB Normalized Reconstruction Error: %e \n',dtu_prob_rec_err);

%compute log spect distance
dtu_prob_logspectdist = mean(sqrt(sum( ...
    reshape( (10*log10((MySTFTabs(:)+eps)./(dtu_prob(:)+eps))).^2 , numrows, numcols)  ...
    ,1)));  
fprintf('      DTU-PROB Final Log Spectral Distance: %e \n',dtu_prob_logspectdist);

% compute BSS EVAL. make sure we define rows as signals, not columns
[dtu_prob_SDR, dtu_prob_SIR, dtu_prob_SAR] = bss_eval_sources(dtu_prob_newsig',origMix');
fprintf('      DTU-PROB SDR: %f \t SIR: %f \t SAR: %f \n',dtu_prob_SDR, dtu_prob_SIR, dtu_prob_SAR);

disp('--- Finished ---')


%% present results in a table for easier incorporation in other docs

T = table(...
[dtu_als_time dtu_alsobs_time dtu_cjlin_time dtu_mm_time dtu_prob_time ]',...
[dtu_als_W_neg dtu_alsobs_W_neg dtu_cjlin_W_neg dtu_mm_W_neg dtu_prob_W_neg ]',...
[dtu_als_H_neg dtu_alsobs_H_neg dtu_cjlin_H_neg dtu_mm_H_neg dtu_prob_H_neg ]',...
[dtu_als_rec_err dtu_alsobs_rec_err dtu_cjlin_rec_err dtu_mm_rec_err dtu_prob_rec_err ]',...
[dtu_als_logspectdist dtu_alsobs_logspectdist dtu_cjlin_logspectdist dtu_mm_logspectdist dtu_prob_logspectdist ]',...
[dtu_als_SDR dtu_alsobs_SDR dtu_cjlin_SDR dtu_mm_SDR dtu_prob_SDR ]',...
'VariableNames', {'Time', 'NegW', 'NegH', 'RecErr', 'LogSpectDist', 'SDR'}, ...
'RowNames', {'DTU ALS','DTU ALSOBS', 'DTU CJLIN', 'DTU MM', 'DTU PROB'}...
);

%% delete unnecessary rows/columns
% like those with negative values
toDelete = T.NegW > 0;
T(toDelete, :) = [];
% now delete NegW and NegH columns
T(:,{'NegW', 'NegH'}) = [];
% display table
T
