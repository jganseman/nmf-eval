% Testing divergences and powers of Beta-NMF - a third experiment. 
% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% october 2017

% For the STFT, we choose the corrected version from CATbox, which is available
% in github.com/jganseman/nsgt_eval . 

% Using the Beta-NMF implementation from G. Grindlay's NMFLIB.
% Phase is ignored in the NMF; the inverse STFT uses the original phase data.
% For separation instead of signal reconstruction, BSS_EVAL 3.0 is used.

% Dataset tested is MAPS piano sounds, (NSynth Test todo?), to create 
% 3-different-note mixtures that are to be separated. Limiting ourselves
% to reasonable MIDI range and fairly equal loudness.

% In this third version of the experiment, we find an average performance
% across 10 tries, using less iterations and a shorter signal

% NOTE: full experiment takes about 1 hour / chord / set of parallel tries.

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
iter = 20; % number of EM iterations
beta = 1;
power = 1;
numtrials = 20;

% define limits of the experiment
lowestpower = 0.5;
powerinc = 0.1; %1.0
highestpower = 2.5;
lowestbeta = 0.0;
highestbeta = 2.5; %2.0;
betainc = 0.1; %1.0;
numpowers = (highestpower - lowestpower) / powerinc +1;
numbetas = (highestbeta - lowestbeta) / betainc +1;
mySDRs = cell([numtrials numpowers numbetas]);
mySIRs = cell([numtrials numpowers numbetas]);
mySARs = cell([numtrials numpowers numbetas]);
myPerm = cell([numtrials numpowers numbetas]);

% pre-initializing other arrays and cell arrays for speed
source = cell([1 nrsources]);
arpegsrc = cell([1 nrsources]);
wavstft = cell([1 nrsources]);
wavstftabs = cell([1 nrsources]);
wavstftphase = cell([1 nrsources]);
reconstnmf = cell([numtrials nrsources]);
nmfrec = cell([numtrials nrsources]);

W0 = cell([1 numtrials]);
H0 = cell([1 numtrials]);

%% begin a set of external loops

for midibase = [60 48 72] % 84 36 96 54 66 42 78 90]
%for thischord= 1:10

% in case of testing, use:
%for midibase = 72
for thischord= 1

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

% create arpeggiated chord plus mixture
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
[nrrows, nrcols] = size(mixstftabs);

%make a set of initial matrices
for trial = 1:numtrials
    W0{trial} = rand(nrrows, nrcomp); %W0 = W0./sum(W0, 1);
    H0{trial} = rand(nrcomp, nrcols); %H0 = H0./sum(H0, 1);
% IMPORTANT: DO NOT MAKE UNIFORM MATRICES. Then the columns can converge
% to the same value if the NMF algo is deterministic. Instead, assure
% linear independence of starting vectors in W0 -> rand
end

% Iterate over power and beta
powerindex=0;

for power = lowestpower:powerinc:highestpower %1
powerindex = powerindex+1;
betaindex=0;

for beta = lowestbeta:betainc:highestbeta %1
betaindex=betaindex+1;   

%print time
fprintf('Starting %d NMF trials, power %f and beta %f. ', numtrials, power, beta);
tic;

% RUN nmf_beta from Grindlay's NMFlab
wavstftabs = origwavstft; % reload original data  
mixstftabs = origmixstft.^power;

% Do x different NMF trials in parallel
parfor trial=1:numtrials

    % perform NMFs
    [nmfW{trial},nmfH{trial}] = nmf_beta(mixstftabs, nrcomp,'niter', iter, 'beta', beta, 'W0', W0{trial}, 'H0', H0{trial});
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

    for i=1:nrcomp
        % Reconstruct sources through Generalized Wiener filtering (DOI 10.1109/ICASSP.2015.7177973)
        reconstnmf{trial, i} = (nmfW{trial}(:,i)*nmfH{trial}(i,:)) ./ (nmfW{trial}*nmfH{trial}); % make mask
        reconstnmf{trial, i} = reconstnmf{trial, i} .* mixstft; % apply to (non-raised) spectrogram
        
        % Inverse STFT. parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
        nmfrec{trial, i} = istft_catbox( reconstnmf{trial, i}, sz/hp, sz, 'smooth' );
        nmfrec{trial, i} = nmfrec{trial, i}(1:size(arpegmix,2));   % cut to size
    end
    
end %parfor NMFs

%print time
timenmf = toc;
fprintf('- Took %f seconds.\n', timenmf);
fprintf('Computing %d BSS_EVAL metrics. ', numtrials);
tic;

% Do x different BSS_EVAL calculations in parallel
% Can not be merged in previous loop for some strange reason
parfor trial=1:numtrials

    [mySDRs{trial, powerindex, betaindex}, mySIRs{trial, powerindex, betaindex}, mySARs{trial, powerindex, betaindex}, myPerm{trial, powerindex, betaindex}] =  bss_eval_sources( cat(1, nmfrec{trial, :}), cell2mat(arpegsrc'));
    %fprintf('Trial %d SDR of sources 1, 2, 3: %f , %f, %f\n', trial, mySDRs{trial, powerindex, betaindex}(1), mySDRs{trial, powerindex, betaindex}(2), mySDRs{trial, powerindex, betaindex}(3));
    %NOTE: BSS_EVAL 3.0 outputs metrics that are already ORDERED according to the TARGET!

end %parfor BSS_EVAL    
timebss = toc;
fprintf('- Took %f seconds.\n', timebss);

end %beta
end %power

%% compute statistics on all trials (manually to make sure we sum correct dimensions)

powerindex=0;
for power = lowestpower:powerinc:highestpower %1
powerindex = powerindex+1;
betaindex=0;
for beta = lowestbeta:betainc:highestbeta %1
betaindex=betaindex+1;    
    % means : concat each trial in column, calc mean over rows
    mySDRmean{powerindex, betaindex} = mean( cat(2, mySDRs{:, powerindex, betaindex}) , 2);
    mySIRmean{powerindex, betaindex} = mean( cat(2, mySIRs{:, powerindex, betaindex}) , 2);
    mySARmean{powerindex, betaindex} = mean( cat(2, mySARs{:, powerindex, betaindex}) , 2);
    % medians
    mySDRmedian{powerindex, betaindex} = median( cat(2, mySDRs{:, powerindex, betaindex}) , 2);
    mySIRmedian{powerindex, betaindex} = median( cat(2, mySIRs{:, powerindex, betaindex}) , 2);
    mySARmedian{powerindex, betaindex} = median( cat(2, mySARs{:, powerindex, betaindex}) , 2);
    % stddevs
    mySDRstddev{powerindex, betaindex} = std( cat(2, mySDRs{:, powerindex, betaindex}) , 0, 2);
    mySIRstddev{powerindex, betaindex} = std( cat(2, mySIRs{:, powerindex, betaindex}) , 0, 2);
    mySARstddev{powerindex, betaindex} = std( cat(2, mySARs{:, powerindex, betaindex}) , 0, 2);
    % max
    mySDRmax{powerindex, betaindex} = max( cat(2, mySDRs{:, powerindex, betaindex}) , 2);
    mySIRmax{powerindex, betaindex} = max( cat(2, mySIRs{:, powerindex, betaindex}) , 2);
    mySARmax{powerindex, betaindex} = max( cat(2, mySARs{:, powerindex, betaindex}) , 2);    
    % min
    mySDRmin{powerindex, betaindex} = min( cat(2, mySDRs{:, powerindex, betaindex}) , 2);
    mySIRmin{powerindex, betaindex} = min( cat(2, mySIRs{:, powerindex, betaindex}) , 2);
    mySARmin{powerindex, betaindex} = min( cat(2, mySARs{:, powerindex, betaindex}) , 2);  
    
end %beta
end %power

%default to mean:
mySDR = mySDRmean;
mySIR = mySIRmean;
mySAR = mySARmean;

%% save to file
disp('Saving to file');
myfilename = ['nmf-beta-test4-chord' int2str(thischord) '-base' int2str(midibase) '.mat']; 
save(myfilename, 'arpegmix', 'arpegsrc', 'chords', 'datadir', 'hp', 'sz', 'iter', 'midibase', ...
    'minlength', 'mixture', 'nrcomp', 'nrcols', 'nrrows', 'nrsources', 'smap', ...
    'mySARmean', 'mySDRmean', 'mySIRmean', 'mySARmedian', 'mySDRmedian', 'mySIRmedian', ...
    'mySARmax', 'mySDRmax', 'mySIRmax', 'mySARmin', 'mySDRmin', 'mySIRmin', ...
    'mySARstddev', 'mySDRstddev', 'mySIRstddev', 'mySAR', 'mySDR', 'mySIR', ...
    'filename');

end %thischord
end %midibase
