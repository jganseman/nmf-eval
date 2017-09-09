% Testing divergences and powers of Beta-NMF - a few experiments 
% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% june 2017

% For the STFT, we choose the corrected version from CATbox, which is available
% in github.com/jganseman/nsgt_eval . This is the fastest STFT among the 4 
% we tested, and also the one with the smallest reconstruction error when
% window-corrected (option 'smooth').

% Using the Beta-NMF implementation from G. Grindlay's NMFLIB
% Phase is ignored in the NMF; the inverse STFT uses the original phase data.
% For separation instead of signal reconstruction, BSS_EVAL 2.1 is used.

% Datasets tested are MAPS piano sounds, and NSynth Test subset, to create 
% 3-different-note mixtures that are to be separated. Limiting ourselves
% to reasonable MIDI range and fairly equal loudness.

% Much of the code in this file is based on an earlier 2014 experiment to 
% run the SISAL algorithm on the MAPS database.

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
beta = 1;
power = 1;

%% begin a set of external loops

for midibase = [60 48 72] % 84 36 96 54 66 42 78 90]
%for thischord= 1:10

% in case of testing, use
%for midibase = 60
for thischord= 1

%% read source files
currentchord = chords{thischord} + midibase;
nrsources = length(currentchord);

minlength = 1000000000;
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
   if length(source{i}) < minlength
       minlength = length(source{i});
   end
end
% Note: minlength is much longer for cases where all sources use sustain=1
% This is the case for: 60-63-68, 60-64-69, 60-64-68
% The inverse (all S=0) is true for: 72-76-79, 72-76-80, 

% cut sources such that they have the same length
mixture=0;
for i=1:nrsources
    source{i} = source{i}(1:minlength);
    mixture = mixture + source{i};
end
mixture = mixture ./ nrsources;

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

parfor i=1:nrsources
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


%% Iterate over power and beta
powerindex=0;

for power = 0.5:0.1:2.5 %1
powerindex = powerindex+1;
betaindex=0;

for beta = 0.0:0.1:2.5 %1
betaindex=betaindex+1;   

fprintf('doing NMF with power %f and beta %f\n', power, beta);

%% RUN nmf_beta from Grindlay's NMFlab

% reload original data    
wavstftabs = origwavstft;
mixstftabs = origmixstft.^power;
%TODO iterate over power and beta

[nmfW,nmfH] = nmf_beta(mixstftabs, nrcomp,'niter', iter, 'beta', beta, 'W0', W0, 'H0', H0);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

%% Reconstruct sources through Generalized Wiener filtering
% See [Liutkus, 2015] DOI 10.1109/ICASSP.2015.7177973 for justification.

parfor i=1:nrcomp
    % first, create all masks  
    reconstnmf{i} = (nmfW(:,i)*nmfH(i,:)) ./ (nmfW*nmfH);
    % apply (raised) masks to original (non-raised) spectrogram
    reconstnmf{i} = reconstnmf{i} .* mixstft;
end

%% TODO sort spectrograms? Maybe best to correlate the reconstructed sources...
% or just pick the highest coefficient of all possible source/reconstr combinations?
% see also https://www.researchgate.net/post/how_to_measure_the_similarity_between_two_signal
% --> maybe correlate the average spectra in frames defined

%% Inverse STFT
disp('Inverse STFT calculation');
% parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
parfor i=1:nrcomp 
   nmfrec{i} = istft_catbox( reconstnmf{i}, sz/hp, sz, 'smooth' );
   nmfrec{i} = nmfrec{i}(1:size(arpegmix,2));   % cut to size
   % find the index of the current source by correlating with part of
   % mixture containing only solo voices
   [ ~ , nmfindex{i} ] = max( nmfrec{i}(1:nrcomp*minlength) .* arpegmix(1:nrcomp*minlength) );
   nmfindex{i} = floor(nmfindex{i}/minlength)+1;
end

% test that nmfindex is a permutation of 1:nrcomp
if ~ismember(cell2mat(nmfindex), perms(1:nrcomp), 'rows')
   fprintf('Could not assign extracted signals to original sources. Skipping.\n');
   % fill up the collected data with empty values
   savednmfrec{powerindex}{betaindex} = cell(size(nmfrec));
   mySDR{powerindex}{betaindex} = [0;0;0];
   mySIR{powerindex}{betaindex} = [0;0;0];
   mySAR{powerindex}{betaindex} = [0;0;0];
   continue;    % break loop, go to next chord/beta/power 
end

%now reorder according to correct indices.
nmfrec(cell2mat(nmfindex)) = nmfrec;
% save this to be able to inspect results later.
savednmfrec{powerindex}{betaindex} = nmfrec;

%% compute the BSS-EVAL SDR metrics
disp('computing BSS-EVAL');

[mySDR{powerindex}{betaindex} mySIR{powerindex}{betaindex} mySAR{powerindex}{betaindex}] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));

fprintf('SDR of sources 1, 2, 3: %f , %f, %f\n', mySDR{powerindex}{betaindex}(1), mySDR{powerindex}{betaindex}(2), mySDR{powerindex}{betaindex}(3));

end %beta
end %power

%% save to file
disp('Saving to file');
myfilename = ['nmf-beta-test-chord' int2str(thischord) '-base' int2str(midibase) '.mat']; 
save(myfilename, 'arpegmix', 'arpegsrc', 'chords', 'datadir', 'H0', 'W0', 'hp', 'sz', 'iter', 'midibase', 'minlength', 'mixture', 'mySAR', 'mySDR', 'mySIR', 'nrcomp', 'nrcols', 'nrrows', 'nrsources', 'smap', 'filename');

    
end %thischord
end %midibase
