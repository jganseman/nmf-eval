% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% november 2017

% This demo tests the convergence in SDR, log spectral distance, or reconstruction error
% of the main different methods to perform NMF: multiplicative updates,
% projected gradient, alternating (nonnegative) least squares, etc.

% For the STFT, we choose the corrected version from CATbox, which is available
% in github.com/jganseman/nsgt_eval . This is the fastest STFT among the 4 
% we tested, and also the one with the smallest reconstruction error when
% window-corrected (option 'smooth').

% Phase is ignored in the NMF; the inverse STFT uses the original phase data.
% For separation instead of signal reconstruction, BSS_EVAL 3.0 is used.
% The separation routines were modified for this purpose and have
% *_verbose attached to their filenames.

%%
clear;
addpath(genpath('../'));     % add subdirectories to path. 
    % Make sure you are in the '/tests' folder, not in the nmf-eval root!
    
%% Starting MAPS experiment
datadir = getMAPS();
if ispc
    datadir = strcat(datadir,'MAPS\ENSTDkAm\ISOL\NO\');
else
    datadir = strcat(datadir,'MAPS/ENSTDkAm/ISOL/NO/');
end
% this is the directory for the Yamaha disklavier with Ambient sound.

fprintf('Starting MAPS BetaNMF experiment, data folder: %s\n', datadir);
    
%% Setup of MAPS filename structure to iterate over files later

%define 3-note chords starting from a root note
chords = { [0 4 7], [0 3 8], [0 5 9], ...       % major chords and inversions
            [0 3 7], [0 4 9], [0 5 8], ...      % minor chords and inversions
            [0 3 6], [0 3 9], [0 6 9], ...      % diminished chords and inversions
            [0 4 8] };                          % augmented chord
    %midibase = 60;      % loop over notes later on.
        
% for 4/5-note chords: take table from MAPS docs (check for duplicates)

% files are called: MAPS_ISOL_NO_x_Sy_Mz_ENSTDkAm.wav, where:
% - x is M, P or F (loudness)
% - y is 0 or 1 (sustain pedal used or not)
% - z is the midi note number, from 21 to 108
% let's limit ourselves here to F notes with no sustain. Ambient recording 
% noise is enough trouble as it is.

filestart = 'MAPS_ISOL_NO_F_S';
filemiddle = '_M';
fileend = '_ENSTDkAm.wav';
%BUT: not all midi nrs present! Some with S0, others with S1 ...
% solution: create a sustain pedal mapping from MIDI note nr to 0/1
smap = zeros(1, 108);        % ugly code, i know, but quick to write and does the trick
              smap(21) = 1; smap(22) = 1; smap(23) = 0; smap(24) = 0; smap(25) = 0; smap(26) = 0; smap(27) = 1; smap(28) = 0; smap(29) = 0; 
smap(30) = 1; smap(31) = 1; smap(32) = 0; smap(33) = 0; smap(34) = 0; smap(35) = 0; smap(36) = 1; smap(37) = 0; smap(38) = 0; smap(39) = 1; 
smap(40) = 0; smap(41) = 1; smap(42) = 0; smap(43) = 1; smap(44) = 1; smap(45) = 0; smap(46) = 1; smap(47) = 1; smap(48) = 1; smap(49) = 0; 
smap(50) = 0; smap(51) = 0; smap(52) = 0; smap(53) = 0; smap(54) = 0; smap(55) = 1; smap(56) = 1; smap(57) = 0; smap(58) = 1; smap(59) = 1; 
smap(60) = 1; smap(61) = 0; smap(62) = 1; smap(63) = 1; smap(64) = 1; smap(65) = 0; smap(66) = 0; smap(67) = 0; smap(68) = 1; smap(69) = 1; 
smap(70) = 0; smap(71) = 1; smap(72) = 0; smap(73) = 0; smap(74) = 1; smap(75) = 1; smap(76) = 0; smap(77) = 1; smap(78) = 1; smap(79) = 0; 
smap(80) = 0; smap(81) = 1; smap(82) = 1; smap(83) = 0; smap(84) = 1; smap(85) = 1; smap(86) = 1; smap(87) = 0; smap(88) = 1; smap(89) = 1; 
smap(90) = 0; smap(91) = 0; smap(92) = 1; smap(93) = 0; smap(94) = 0; smap(95) = 1; smap(96) = 1; smap(97) = 0; smap(98) = 0; smap(99) = 0; 
smap(100) = 0; smap(31) = 0; smap(102) = 0; smap(103) = 0; smap(104) = 0; smap(105) = 1; smap(106) = 0; smap(107) = 0; smap(108) = 1; 

%% Setup of STFT and NMF parameters
nrsources=3;
nrcomp = 3;
sz = 1024; % stft frame size
hp = sz / 4; % hop size
iter = 50; % number of EM iterations
power = 1;

midibase = 60;
thischord= 1;

%% read source files
currentchord = chords{thischord} + midibase;
nrsources = length(currentchord);

minlength = 128000/8;     % 0.75 seconds is more than enough
fprintf('Processing files: ')
for i=1:nrsources
   filename = [ filestart num2str(smap(currentchord(i))) filemiddle num2str(currentchord(i)) fileend ]; 
   fprintf('%s ; ', filename);
   source{i} = audioread([ datadir filename ]);
   % in case s1 is used, the added pedaling makes that the file is 0.2
   % seconds longer in the beginning. Cut this extra time off.
   if smap(currentchord(i))     %equals 1 is implicit
       source{i} = source{i}(0.2*44100:end);
   end
end
% Note: minlength is much longer for cases where all sources use sustain=1
% This is the case for: 60-63-68, 60-64-69, 60-64-68
% The inverse (all S=0) is true for: 72-76-79, 72-76-80, 

% cut sources such that they have the same length
mixture=0;
for i=1:nrsources
    source{i} = source{i}(22050+1:22050+minlength);
    mixture = mixture + source{i};
end
mixture = mixture ./ nrsources;
%NOTE: there was still a +/- 0.5s silence in beginning of each source!

% now we've got our sources in the source{} cell, and a mixture.
% we can concatenate the sources to our mixture, to have an arpeggiated
% chord plus the mixture. for sdr comparison later, make new sources as well
arpegmix = [ cell2mat(source) mixture ];
silence=zeros(1,minlength);
%the next lines are only for 3-note chords. adapt for 4 note!
arpegsrc{1} = [ source{1} silence silence source{1}./nrsources ];
arpegsrc{2} = [ silence source{2} silence source{2}./nrsources ];
arpegsrc{3} = [ silence silence source{3} source{3}./nrsources ];

fprintf('\nLength of mixture: %2.3f seconds.\n', minlength*4/44100);
% to play, use:
% soundsc(arpegmix, 44100)

%% compute STFTs

disp('--- Computing corrected CATbox STFTs with periodic Hann window ---');  
% parameters: signal, window, overlap, fftsize

for i=1:nrsources
    % spectrograms, magnitudes and phases
    wavstft{i} = stft_catbox(arpegsrc{i}, hann(sz, 'periodic'), sz-hp, sz);
    wavstftabs{i} = max(abs(wavstft{i}), eps);
    wavstftphase{i} = wavstft{i} ./ wavstftabs{i};
end
    mixstft = stft_catbox(arpegmix, hann(sz, 'periodic'), sz-hp, sz);
    mixstftabs = max(abs(mixstft), eps);
    mixstftphase = mixstft ./ mixstftabs;
    
% make original copy of stfts, to reload before each computation
origwavstft = wavstftabs;
origmixstft = mixstftabs;

% create starting matrices;
[nrrows, nrcols] = size(mixstftabs);
% IMPORTANT: DO NOT MAKE UNIFORM MATRICES. Then the columns can converge
% to the same value if the NMF algo is deterministic. Instead, assure
% linear independence of starting vectors in W0 -> rand
W0 = rand(nrrows, nrcomp); %W0 = W0./sum(W0, 1);
H0 = rand(nrcomp, nrcols); %H0 = H0./sum(H0, 1);


%% RUN nmf_beta from Grindlay's NMFlab
disp('Testing NMFLAB beta (mult. upd., KL-div)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
beta = 1;
[beta1W, beta1H, beta1_errs, beta1_iterW, beta1_iterH, beta1_lsd] = nmf_beta_verbose(mixstftabs, nrcomp,'niter', iter, 'beta', beta, 'W0', W0, 'H0', H0);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

beta1_spec_err = zeros(iter,1);
beta1_fullreconst = zeros(iter, size(arpegmix,2));
beta1_signal_err = zeros(iter,1);
beta1_sdr = zeros(nrcomp, iter);
beta1_sir = zeros(nrcomp, iter);
beta1_sar = zeros(nrcomp, iter);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    beta1_spec_err(j) = norm(mixstftabs- (beta1_iterW(:,:,j)*beta1_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    beta1_temp = istft_catbox( (beta1_iterW(:,:,j)*beta1_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    beta1_fullreconst(j,:) = beta1_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    beta1_signal_err(j) = norm(arpegmix - beta1_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources 
    for i=1:nrcomp
        % first, create all masks  
        reconstbeta1{i} = (beta1_iterW(:,i,j)*beta1_iterH(i,:,j)) ./ (beta1_iterW(:,:,j)*beta1_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstbeta1{i} = reconstbeta1{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       beta1rec{i} = istft_catbox( reconstbeta1{i}, sz/hp, sz, 'smooth' );
       beta1rec{i} = beta1rec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [beta1_sdr(:,j), beta1_sir(:,j), beta1_sar(:,j)] =  bss_eval_sources(cell2mat(beta1rec'), cell2mat(arpegsrc')); 
end


%% RUN nmf_beta from Grindlay's NMFlab
disp('Testing NMFLAB beta (mult. upd., IS-div)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
beta = 0;
[beta0W, beta0H, beta0_errs, beta0_iterW, beta0_iterH, beta0_lsd] = nmf_beta_verbose(mixstftabs, nrcomp,'niter', iter, 'beta', beta, 'W0', W0, 'H0', H0);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

beta0_spec_err = zeros(iter,1);
beta0_fullreconst = zeros(iter, size(arpegmix,2));
beta0_signal_err = zeros(iter,1);
beta0_sdr = zeros(nrcomp, iter);
beta0_sir = zeros(nrcomp, iter);
beta0_sar = zeros(nrcomp, iter);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    beta0_spec_err(j) = norm(mixstftabs- (beta0_iterW(:,:,j)*beta0_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    beta0_temp = istft_catbox( (beta0_iterW(:,:,j)*beta0_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    beta0_fullreconst(j,:) = beta0_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    beta0_signal_err(j) = norm(arpegmix - beta0_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources 
    for i=1:nrcomp
        % first, create all masks  
        reconstbeta0{i} = (beta0_iterW(:,i,j)*beta0_iterH(i,:,j)) ./ (beta0_iterW(:,:,j)*beta0_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstbeta0{i} = reconstbeta0{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       beta0rec{i} = istft_catbox( reconstbeta0{i}, sz/hp, sz, 'smooth' );
       beta0rec{i} = beta0rec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [beta0_sdr(:,j), beta0_sir(:,j), beta0_sar(:,j)] =  bss_eval_sources(cell2mat(beta0rec'), cell2mat(arpegsrc')); 
end


%% RUN nmf_beta from Grindlay's NMFlab
disp('Testing NMFLAB beta (mult. upd., EUC-div)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
beta = 2;
[beta2W, beta2H, beta2_errs, beta2_iterW, beta2_iterH, beta2_lsd] = nmf_beta_verbose(mixstftabs, nrcomp,'niter', iter, 'beta', beta, 'W0', W0, 'H0', H0);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

beta2_spec_err = zeros(iter,1);
beta2_fullreconst = zeros(iter, size(arpegmix,2));
beta2_signal_err = zeros(iter,1);
beta2_sdr = zeros(nrcomp, iter);
beta2_sir = zeros(nrcomp, iter);
beta2_sar = zeros(nrcomp, iter);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    beta2_spec_err(j) = norm(mixstftabs- (beta2_iterW(:,:,j)*beta2_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    beta2_temp = istft_catbox( (beta2_iterW(:,:,j)*beta2_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    beta2_fullreconst(j,:) = beta2_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    beta2_signal_err(j) = norm(arpegmix - beta2_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources 
    for i=1:nrcomp
        % first, create all masks  
        reconstbeta2{i} = (beta2_iterW(:,i,j)*beta2_iterH(i,:,j)) ./ (beta2_iterW(:,:,j)*beta2_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstbeta2{i} = reconstbeta2{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       beta2rec{i} = istft_catbox( reconstbeta2{i}, sz/hp, sz, 'smooth' );
       beta2rec{i} = beta2rec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [beta2_sdr(:,j), beta2_sir(:,j), beta2_sar(:,j)] =  bss_eval_sources(cell2mat(beta2rec'), cell2mat(arpegsrc')); 
end




%% RUN nmf_als from Hansen's NMF_DTU
disp('Testing NMF_DTU ALS (regular least squares)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
[alsW,alsH,als_errs, als_iterW, als_iterH, als_lsd] = nmf_als_verbose(mixstftabs, nrcomp,'Nitsmax', iter, 'W0', W0, 'H0', H0);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

als_spec_err = zeros(iter,1);
als_sdr = zeros(nrcomp, iter);
als_sir = zeros(nrcomp, iter);
als_sar = zeros(nrcomp, iter);
als_fullreconst = zeros(iter, size(arpegmix,2));
als_signal_err = zeros(iter,1);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    als_spec_err(j) = norm(mixstftabs- (als_iterW(:,:,j)*als_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    als_temp = istft_catbox( (als_iterW(:,:,j)*als_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    als_fullreconst(j,:) = als_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    als_signal_err(j) = norm(arpegmix - als_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources
    for i=1:nrcomp
        % first, create all masks  
        reconstals{i} = (als_iterW(:,i,j)*als_iterH(i,:,j)) ./ (als_iterW(:,:,j)*als_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstals{i} = reconstals{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       alsrec{i} = istft_catbox( reconstals{i}, sz/hp, sz, 'smooth' );
       alsrec{i} = alsrec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [als_sdr(:,j), als_sir(:,j), als_sar(:,j)] =  bss_eval_sources(cell2mat(alsrec'), cell2mat(arpegsrc')); 
end


%% RUN nmf_cjlin from Hansen's NMF_DTU
disp('Testing NMF_DTU CJLIN (gradient descent)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
[cjlinW,cjlinH,cjlin_errs, cjlin_iterW, cjlin_iterH, cjlin_lsd] = nmf_cjlin_verbose(mixstftabs,'Winit', W0, 'Hinit', H0, 'tol', 10e-5, 'maxiter', iter);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

cjlin_spec_err = zeros(iter,1);
cjlin_sdr = zeros(nrcomp, iter);
cjlin_sir = zeros(nrcomp, iter);
cjlin_sar = zeros(nrcomp, iter);
cjlin_fullreconst = zeros(iter, size(arpegmix,2));
cjlin_signal_err = zeros(iter,1);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    cjlin_spec_err(j) = norm(mixstftabs- (cjlin_iterW(:,:,j)*cjlin_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    cjlin_temp = istft_catbox( (cjlin_iterW(:,:,j)*cjlin_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    cjlin_fullreconst(j,:) = cjlin_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    cjlin_signal_err(j) = norm(arpegmix - cjlin_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources
    for i=1:nrcomp
        % first, create all masks  
        reconstcjlin{i} = (cjlin_iterW(:,i,j)*cjlin_iterH(i,:,j)) ./ (cjlin_iterW(:,:,j)*cjlin_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstcjlin{i} = reconstcjlin{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       cjlinrec{i} = istft_catbox( reconstcjlin{i}, sz/hp, sz, 'smooth' );
       cjlinrec{i} = cjlinrec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [cjlin_sdr(:,j), cjlin_sir(:,j), cjlin_sar(:,j)] =  bss_eval_sources(cell2mat(cjlinrec'), cell2mat(arpegsrc')); 
end



%% RUN blockpivot from NMF_APD
disp('Testing NMF_APD BlockPivot (Nonnegative Least Squares)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
opts.W0 = W0; opts.H0 = H0; 
opts.maxIter = iter;
opts.verbose = 0;
[bpW,bpH,bp_errs, bp_iterW, bp_iterH, bp_lsd] = nmf_bp_verbose(mixstftabs, opts, 'NNLS_SOLVER', 'bp', 'MAX_ITER', iter, 'TOL', eps, 'VERBOSE', 0);

bp_spec_err = zeros(iter,1);
bp_sdr = zeros(nrcomp, iter);
bp_sir = zeros(nrcomp, iter);
bp_sar = zeros(nrcomp, iter);
bp_fullreconst = zeros(iter, size(arpegmix,2));
bp_signal_err = zeros(iter,1);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    bp_spec_err(j) = norm(mixstftabs- (bp_iterW(:,:,j)*bp_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    bp_temp = istft_catbox( (bp_iterW(:,:,j)*bp_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    bp_fullreconst(j,:) = bp_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    bp_signal_err(j) = norm(arpegmix - bp_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources
    for i=1:nrcomp
        % first, create all masks  
        reconstbp{i} = (bp_iterW(:,i,j)*bp_iterH(i,:,j)) ./ (bp_iterW(:,:,j)*bp_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstbp{i} = reconstbp{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       bprec{i} = istft_catbox( reconstbp{i}, sz/hp, sz, 'smooth' );
       bprec{i} = bprec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [bp_sdr(:,j), bp_sir(:,j), bp_sar(:,j)] =  bss_eval_sources(cell2mat(bprec'), cell2mat(arpegsrc')); 
end



%% RUN activeSet from NMF_APD
% Note: seems to compute same solution as BlockPivot, hence don't compare.

disp('Testing NMF_APD ActiveSet (Nonnegative Least Squares)');
% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;

% compute NMF and save history of variables
opts.W0 = W0; opts.H0 = H0; 
opts.maxIter = iter;
opts.verbose = 0;
[asW,asH,as_errs, as_iterW, as_iterH, as_lsd] = nmf_bp_verbose(mixstftabs, opts, 'NNLS_SOLVER', 'as', 'MAX_ITER', iter, 'TOL', eps, 'VERBOSE', 0);

as_spec_err = zeros(iter,1);
as_sdr = zeros(nrcomp, iter);
as_sir = zeros(nrcomp, iter);
as_sar = zeros(nrcomp, iter);
as_fullreconst = zeros(iter, size(arpegmix,2));
as_signal_err = zeros(iter,1);
% now go over all intermediate iterations
for j=1:iter
    fprintf('Reconstructing/evaluating iteration %d\n', j); 
    % Spectral reconstruction error
    as_spec_err(j) = norm(mixstftabs- (as_iterW(:,:,j)*as_iterH(:,:,j)) )/norm(mixstftabs);
    % Reconstruct full signal and cut to size
    as_temp = istft_catbox( (as_iterW(:,:,j)*as_iterH(:,:,j)).*mixstftphase , sz/hp, sz, 'smooth' );
    as_fullreconst(j,:) = as_temp(1:size(arpegmix,2));
    % Total signal reconstruction error
    as_signal_err(j) = norm(arpegmix - as_fullreconst(j,:) )/norm(arpegmix);
    % Reconstruct sources
    for i=1:nrcomp
        % first, create all masks  
        reconstas{i} = (as_iterW(:,i,j)*as_iterH(i,:,j)) ./ (as_iterW(:,:,j)*as_iterH(:,:,j));
        % apply (raised) masks to original (non-raised) spectrogram
        reconstas{i} = reconstas{i} .* mixstft;
    end
    % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
    for i=1:nrcomp 
       asrec{i} = istft_catbox( reconstas{i}, sz/hp, sz, 'smooth' );
       asrec{i} = asrec{i}(1:size(arpegmix,2));   % cut to size
       % see plotBetaNMF.m for code to eventually order the sources.
    end
    % compute the BSS-EVAL SDR metrics
    [as_sdr(:,j), as_sir(:,j), as_sar(:,j)] =  bss_eval_sources(cell2mat(asrec'), cell2mat(arpegsrc')); 
end

    
%% plot SDR of first source

PRINTTOFILE = 1;


% figure; 
% plot(1:1:iter, beta1_sdr(1,:), 'blue'); hold on;
% plot(1:1:iter, als_sdr(1,:), 'red'); hold on;
% plot(1:1:iter, cjlin_sdr(1,:), 'green'); hold on;
% plot(1:1:iter, bp_sdr(1,:), 'black'); 
% title('SDR of source 1');
% xlabel('iteration');
% ylabel('dB');
% grid on;
% legend(...
%     {'NMF-Beta-KL',...    %NMF MU KL:  
%     'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},...   
%     'Interpreter', 'latex', 'Location', 'southeast');

% figure; 
% plot(1:1:iter, beta1_sdr(2,:), 'blue'); hold on;
% plot(1:1:iter, als_sdr(2,:), 'red'); hold on;
% plot(1:1:iter, cjlin_sdr(2,:), 'green'); hold on;
% plot(1:1:iter, bp_sdr(2,:), 'black'); 
% title('SDR of source 2');
% xlabel('iteration');
% ylabel('dB');
% grid on;
% legend(...
%     {'NMF-Beta-KL',...    %NMF MU KL:  
%     'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},...   
%     'Interpreter', 'latex', 'Location', 'southeast');

% figure;
% plot(1:1:iter, beta1_sdr(3,:), 'blue'); hold on;
% plot(1:1:iter, als_sdr(3,:), 'red'); hold on;
% plot(1:1:iter, cjlin_sdr(3,:), 'green'); hold on;
% plot(1:1:iter, bp_sdr(3,:), 'black'); 
% title('SDR of source 3');
% xlabel('iteration');
% ylabel('dB');
% grid on;
% legend(...
%     {'NMF-Beta-KL',...    %NMF MU KL:  
%     'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},...   
%     'Interpreter', 'latex', 'Location', 'southeast');

figure; 
plot(1:1:iter, sum(beta1_sdr, 1)./3, '-b'); hold on;
plot(1:1:iter, sum(beta0_sdr, 1)./3, '--b'); hold on;
plot(1:1:iter, sum(beta2_sdr, 1)./3, ':b'); hold on;
plot(1:1:iter, sum(als_sdr, 1)./3, 'red'); hold on;
plot(1:1:iter, sum(cjlin_sdr, 1)./3, 'green'); hold on;
plot(1:1:iter, sum(bp_sdr, 1)./3, '-k'); %hold on; %black
%plot(1:1:iter, sum(as_sdr, 1)./3, '--k');
%title('Average SDR');
xlabel('iteration');
ylabel('dB');
grid on;
%legend(...
%    {'NMF-Beta-KL','NMF-Beta-IS','NMF-Beta-EUC',...   
%    'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},... %}, 'NMF-ActSet'},...   
%    'Interpreter', 'latex', 'Location', 'southeastoutside');

% if wanted, print to file
if PRINTTOFILE
    filename = '../../../thesis/images/chapter5/evolution-sdr.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

% make a copy with only an excerpt
xlim([0 25]);
ylim([16 18]);
if PRINTTOFILE
    filename = '../../../thesis/images/chapter5/evolution-sdr-detail.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end


% for fun: LOG SDR (might reduce effect of scale difference between sources?)
% figure; 
% plot(1:1:iter, sum(log10(beta1_sdr), 1)./3, 'blue'); hold on;
% plot(1:1:iter, sum(log10(als_sdr), 1)./3, 'red'); hold on;
% plot(1:1:iter, sum(log10(cjlin_sdr), 1)./3, 'green'); hold on;
% plot(1:1:iter, sum(log10(bp_sdr), 1)./3, 'black'); 
% title('Average log(SDR)');
% xlabel('iteration');
% ylabel('log(dB)');
% grid on;
% legend(...
%     {'NMF-Beta-KL',...    %NMF MU KL:  
%     'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},...   
%     'Interpreter', 'latex', 'Location', 'southeast');


% plot spectral reconstruction error
% Note: very similar to signal reconstruction error.

% figure; 
% plot(1:1:iter, beta1_spec_err, '-b'); hold on;
% plot(1:1:iter, beta0_spec_err, '--b'); hold on;
% plot(1:1:iter, beta2_spec_err, ':b'); hold on;
% plot(1:1:iter, als_spec_err, 'red'); hold on;
% plot(1:1:iter, cjlin_spec_err, 'green'); hold on;
% plot(1:1:iter, bp_spec_err, '-k'); hold on; %black
% plot(1:1:iter, as_spec_err, '--k'); 
% title('Spectral Reconstruction Error');
% xlabel('iteration');
% ylabel('total error');
% grid on;
% legend(...
%     {'NMF-Beta-KL','NMF-Beta-IS','NMF-Beta-EUC',...   
%     'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP', 'NMF-ActSet'},...   
%     'Interpreter', 'latex', 'Location', 'southeast');

% plot signal reconstruction error
figure; 
plot(1:1:iter, beta1_signal_err, '-b'); hold on;
plot(1:1:iter, beta0_signal_err, '--b'); hold on;
plot(1:1:iter, beta2_signal_err, ':b'); hold on;
plot(1:1:iter, als_signal_err, 'red'); hold on;
plot(1:1:iter, cjlin_signal_err, 'green'); hold on;
plot(1:1:iter, bp_signal_err, '-k'); %hold on; %black
%plot(1:1:iter, as_signal_err, '--k');  
%title('Signal Reconstruction Error');
xlabel('iteration');
ylabel('total error');
grid on;
%legend(...
%    {'NMF-Beta-KL','NMF-Beta-IS','NMF-Beta-EUC',...   
%    'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},... %, 'NMF-ActSet'},...   
%    'Interpreter', 'latex', 'Location', 'northeastoutside');

% if wanted, print to file
if PRINTTOFILE
    filename = '../../../thesis/images/chapter5/evolution-recerr.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%make a copy with only an excerpt
xlim([0 25]);
ylim([0.1 0.4]);
if PRINTTOFILE
    filename = '../../../thesis/images/chapter5/evolution-recerr-detail.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

% plot log spectral distance
figure; 
plot(1:1:iter, beta1_lsd, '-b'); hold on;
plot(1:1:iter, beta0_lsd, '--b'); hold on;
plot(1:1:iter, beta2_lsd, ':b'); hold on;
plot(1:1:iter, als_lsd, 'red'); hold on;
plot(1:1:iter, cjlin_lsd, 'green'); hold on;
plot(1:1:iter, bp_lsd, '-k'); %hold on; %black
%plot(1:1:iter, as_lsd, '--k');  
%title('Mean Log Spectral Distance');
xlabel('iteration');
ylabel('mean total D_{LS}');
grid on;
legend(...
    {'NMF-Beta-KL','NMF-Beta-IS','NMF-Beta-EUC',...   
    'NMF-ALS', 'NMF-ProjG', 'NMF-BlockP'},... %, 'NMF-ActSet'},...   
    'Interpreter', 'latex', 'Location', 'northeastoutside');

% if wanted, print to file
ylim([50 200])
if PRINTTOFILE
    filename = '../../../thesis/images/chapter5/evolution-lsd.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end


% TODO plot KL-divergence? (as recorded in beta_errs)?

