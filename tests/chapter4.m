% Chapter 4: creating the images for Chapter 4 of my PhD Thesis.
% (c) Joachim Ganseman, september 2017

% NOTE: run this file while inside the 'tests' folder in Matlab for correct path

% NOTE: AFTER RUNNING, EXTERNALIZE THE DATA TABLE READING in 2 steps:
% 1. add the following code before subimporting the created tex file:
%     \pgfplotstableread{./images/chapter2/audiosignalexcerpt-1.tsv}\loadedtable
% 2. in the tex file itself, change the following:
%     table[row sep=crcr,format=inline]{\loadedtable}; 

clear;

%%%% PART ONE: illustrate audio-to-score alignment %%%%
% This is based on [Ganseman 2010 ICMC]

PRINTTOFILE=0;

%%
addpath(genpath('./'));     % add subdirectories to path
sz = 2048; % stft frame size
hp = sz / 4; % hop size        Ellis vocoder: 75% overlap needed for reconstruction
emphasis = 1; % spectrogram emphasis value
    %NOTE: emphasis value is used to increase the energy in higher frequency
    %components. This tends to impose more attention to the harmony, since
    %harmonics are composed of higher frequencies but usually lower in energy


% Read file
datadir = getDataDirectory();       % directory with example files
[d1,samplerate] = wavread([datadir 'air-mono-excerpt.wav']);
[d2,samplerate2] = wavread([datadir 'bach-air-score-excerpt.wav']);
% these files are mono

%use CATBOX stft
D1 = stft_catbox(d1, hann(sz, 'periodic'), sz-hp, sz); % 2048pt dft with 512 hop size
D2 = stft_catbox(d2, hann(sz, 'periodic'), sz-hp, sz);

%emphasize higher harmonics
D1 = D1 .* repmat( linspace( 1, emphasis, rows( D1 ))', 1, cols( D1 ));
D2 = D2 .* repmat( linspace( 1, emphasis, rows( D2 ))', 1, cols( D2 ));

%take absolute value
D1abs = abs(D1);      
D2abs = abs(D2);

% draw the results. Use subplots
figure

subplot(2,1,1)
imagesc(20*log(D1abs));
colormap(jet);
axis xy;
xlabel('time')
ylabel('frequency bin')
%title('Real recording')
%put titles in latex inclusion

subplot(2,1,2)
imagesc(20*log(D2abs));
colormap(jet);
axis xy;
xlabel('time')
ylabel('frequency bin')
%title('Synthesized score')

%save this plot if you want
if PRINTTOFILE
    filename='../../../thesis/images/chapter4/bwv1068spect.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% Perform DTW between these audiofiles
% Code based on http://www.ee.columbia.edu/ln/rosa/matlab/dtw/
% make sure MEX file is compiled in helperfuncs!

% Construct the 'local match' scores matrix as the cosine distance 
% between the STFT magnitudes
SM = simmx(D1abs,D2abs);      %need absolute spectrum for this
% Look at it:
figure
imagesc(SM),
colormap(1-gray), 
axis xy;
xlabel('Score time')
ylabel('Recording time')
%title('Similarity Matrix and DTW Path')
% leave title to figure caption in latex
% need a new figure for this, since each figure window can only have 1 colormap. 

% Use dynamic programming to find the lowest-cost path between the 
% opposite corners of the cost matrix
% Note that we use 1-SM because dp will find the *lowest* total cost
[p,q,C] = dpfast(1-SM);
% Overlay the path on the local similarity matrix
hold on; plot(q,p,'r'); hold off
% Path visibly follows the dark stripe

%save this plot if you want
if PRINTTOFILE
    filename = '../../../thesis/images/chapter4/dtwexample.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end


%% Next: check on different result of BSS_EVAL versions
%setup NMF and STFT parameters
nrsources=2;
nrcomp = 2;
sz = 1024; % stft frame size
hp = sz / 4; % hop size
iter = 50; % number of EM iterations
beta = 1;
power = 1;

% use a basic separation example here
datadir = getMAPS();
if ispc
    datadir = strcat(datadir,'MAPS\ENSTDkAm\ISOL\NO\');
else
    datadir = strcat(datadir,'MAPS/ENSTDkAm/ISOL/NO/');
end
% this is the directory for the Yamaha disklavier with Ambient sound.

% read sound files of 2 notes without pedal. Taking MIDI nrs 61 and 67 or 70
source{1} = audioread([ datadir 'MAPS_ISOL_NO_F_S0_M61_ENSTDkAm.wav' ]);
source{2} = audioread([ datadir 'MAPS_ISOL_NO_F_S0_M70_ENSTDkAm.wav' ]);
% not equal in size, contains lots of silence. Trimming to a multiple of 256
minlength = 128000;
%create mixture
mixture=0;
for i=1:nrsources
    source{i} = source{i}(1:minlength);
    mixture = mixture + source{i};
end
mixture = mixture ./ nrsources;

% make mixture of arpeggiated chord followed by full chord
arpegmix = [ cell2mat(source) mixture ];
silence=zeros(1,minlength);
% the ground truth of that arpegmix variable is created below 
% the next lines are only for 2-note chords. adapt for 4 note!
arpegsrc{1} = [ source{1} silence source{1}./nrsources ];
arpegsrc{2} = [ silence source{2} source{2}./nrsources ];

totlength = minlength*3;
fprintf('\nLength of mixture: %2.3f seconds.\n', totlength/44100);
% to play, use:
% soundsc(arpegmix, 44100)

%% compute stfts and run separation
for i=1:nrsources
    % spectrograms, magnitudes and phases
    wavstft{i} = stft_catbox(arpegsrc{i}, hann(sz, 'periodic'), sz-hp, sz);
    wavstftabs{i} = max(abs(wavstft{i}), eps);
    wavstftphase{i} = wavstft{i} ./ wavstftabs{i};
end
mixstft = stft_catbox(arpegmix, hann(sz, 'periodic'), sz-hp, sz);
mixstftabs = max(abs(mixstft), eps);
mixstftphase = mixstft ./ mixstftabs;
    
[nrrows, nrcols] = size(mixstftabs);
W0 = rand(nrrows, nrcomp); %W0 = W0./sum(W0, 1);
H0 = rand(nrcomp, nrcols); %H0 = H0./sum(H0, 1);

% use grindlays nmf
[nmfW,nmfH] = nmf_beta(mixstftabs, nrcomp,'niter', iter, 'beta', beta, 'W0', W0, 'H0', H0);
    % other options: 'thresh', [], 'norm_w', 1, 'norm_h', 0, 'verb', 1, 'myeps', 1e-20, 'W', [], 'H', []);

% reconstruct sources
for i=1:nrcomp
    reconstnmf{i} = (nmfW(:,i)*nmfH(i,:)) ./ (nmfW*nmfH);
    % apply as mask to original spectrogram
    reconstnmf{i} = reconstnmf{i} .* mixstft;
end

% inverse stft
for i=1:nrcomp 
   nmfrec{i} = istft_catbox( reconstnmf{i}, sz/hp, sz, 'smooth' );
   nmfrec{i} = nmfrec{i}(1:size(arpegmix,2));   % cut to size
   % find the index of the current source by correlating with part of mixture containing only solo voices
   [ ~ , nmfindex{i} ] = max( nmfrec{i}(1:nrcomp*minlength) .* arpegmix(1:nrcomp*minlength) );
   nmfindex{i} = floor(nmfindex{i}/minlength)+1;
end

% check if each source is assigned
if ~ismember(cell2mat(nmfindex), perms(1:nrcomp), 'rows')
   fprintf('Could not assign extracted signals to original sources. Skipping.\n');
end

% now reorder according to correct indices.
nmfrec(cell2mat(nmfindex)) = nmfrec;
% To play, use:
% soundsc(nmfrec{1}, 44100)
% soundsc(nmfrec{2}, 44100)


%% compute BSS_EVAL 2.1 sources metrics
[s_target1, e_interf1, e_artif1] = bss_decomp_gain(nmfrec{1}, 1, cell2mat(arpegsrc'));
[sdr21g(1), sir21g(1), sar21g(1)] =  bss_crit(s_target1, e_interf1, e_artif1);
[s_target2, e_interf2, e_artif2] = bss_decomp_gain(nmfrec{2}, 2, cell2mat(arpegsrc'));
[sdr21g(2), sir21g(2), sar21g(2)] =  bss_crit(s_target2, e_interf2, e_artif2);
% note: decomp_gain equals decomp_filt!

% For bss_decomp_tvgain and bss_decomp_tvfilt, a filter shape and size
% needs to be defined as parameter to these functions

fprintf('SDR21g of sources 1, 2: %f , %f\n', sdr21g(1), sdr21g(2));
fprintf('SIR21g of sources 1, 2: %f , %f\n', sir21g(1), sir21g(2));
fprintf('SAR21g of sources 1, 2: %f , %f\n\n', sar21g(1), sar21g(2));

%% compute BSS_EVAL 3.0 sources metrics
[sdr30s, sir30s, sar30s, perm30s] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('SDR30s of sources 1, 2: %f , %f\n', sdr30s(1), sdr30s(2));
fprintf('SIR30s of sources 1, 2: %f , %f\n', sir30s(1), sir30s(2));
fprintf('SAR30s of sources 1, 2: %f , %f\n\n', sar30s(1), sar30s(2));

%% compute BSS_EVAL 3.0 images metrics
[sdr30i,isr30i,sir30i,sar30i,perm30i]=bss_eval_images(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('SDR30i of sources 1, 2: %f , %f\n', sdr30i(1), sdr30i(2));
fprintf('SIR30i of sources 1, 2: %f , %f\n', sir30i(1), sir30i(2));
fprintf('SAR30i of sources 1, 2: %f , %f\n', sar30i(1), sar30i(2));
fprintf('ISR30i of sources 1, 2: %f , %f\n\n', isr30i(1), isr30i(2));

%% compute PEASS 2.0 objective metrics
% done by repurposing the PEASS handler file from 2012 Fritsch-Ganseman ICML paper
options.segmentationFactor = 1;
% Segmentation factor of x splits the file in x pieces. 1 is default
options.resultdir = './'; % where to write the temp results. must end in /
options.sr = 44100;
[ sdrP2, sirP2, sarP2, isrP2, qGlobal, qTarget, qInterf, qArtif, OPS, TPS, IPS, APS ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

% in PEASS toolkit itself (1). Test the difference!
fprintf('sdrP2 of sources 1, 2: %f , %f\n', sdrP2(1), sdrP2(2));
fprintf('sirP2 of sources 1, 2: %f , %f\n', sirP2(1), sirP2(2));
fprintf('sarP2 of sources 1, 2: %f , %f\n', sarP2(1), sarP2(2));
fprintf('isrP2 of sources 1, 2: %f , %f\n\n', isrP2(1), isrP2(2));

%% table output for latex
fprintf('$\\hat{s_1}$ SDR: & $%f$ & $%f$ & $%f$ & $\\hat{s_2}$ SDR: & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sdr21g(1), sdr30s(1), sdr30i(1), sdr21g(2), sdr30s(2), sdr30i(2));
fprintf('$\\hat{s_1}$ SIR: & $%f$ & $%f$ & $%f$ & $\\hat{s_2}$ SIR: & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sir21g(1), sir30s(1), sir30i(1), sir21g(2), sir30s(2), sir30i(2));
fprintf('$\\hat{s_1}$ SAR: & $%f$ & $%f$ & $%f$ & $\\hat{s_2}$ SAR: & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sar21g(1), sar30s(1), sar30i(1), sar21g(2), sar30s(2), sar30i(2));

fprintf('$\\hat{s_1}$ SDR: & $%2.2f$ & $%2.2f$ & $%2.2f$ & $\\hat{s_2}$ SDR: & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr21g(1), sdr30s(1), sdr30i(1), sdr21g(2), sdr30s(2), sdr30i(2));
fprintf('$\\hat{s_1}$ SIR: & $%2.2f$ & $%2.2f$ & $%2.2f$ & $\\hat{s_2}$ SIR: & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sir21g(1), sir30s(1), sir30i(1), sir21g(2), sir30s(2), sir30i(2));
fprintf('$\\hat{s_1}$ SAR: & $%2.2f$ & $%2.2f$ & $%2.2f$ & $\\hat{s_2}$ SAR: & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sar21g(1), sar30s(1), sar30i(1), sar21g(2), sar30s(2), sar30i(2));

fprintf('\n');
%% more table output for latex
fprintf('BSS\\_EVAL 2.1 : & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sdr21g(1), sir21g(1), sar21g(1), sdr21g(2), sir21g(2), sar21g(2));
fprintf('BSS\\_EVAL 3.0\\_s : & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sdr30s(1), sir30s(1), sar30s(1), sdr30s(2), sir30s(2), sar30s(2));
fprintf('BSS\\_EVAL 3.0\\_i : & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sdr30i(1), sir30i(1), sar30i(1), sdr30i(2), sir30i(2), sar30i(2));
fprintf('PEASS 2.0 : & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ & $%f$ \\\\\n\\hline\n', sdrP2(1), sirP2(1), sarP2(1), sdrP2(2), sirP2(2), sarP2(2));

fprintf('BSS\\_EVAL 2.1 : & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr21g(1), sir21g(1), sar21g(1), sdr21g(2), sir21g(2), sar21g(2));
fprintf('BSS\\_EVAL 3.0\\_s : & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30s(1), sir30s(1), sar30s(1), sdr30s(2), sir30s(2), sar30s(2));
fprintf('BSS\\_EVAL 3.0\\_i : & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30i(1), sir30i(1), sar30i(1), sdr30i(2), sir30i(2), sar30i(2));
fprintf('PEASS 2.0 : & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2(1), sirP2(1), sarP2(1), sdrP2(2), sirP2(2), sarP2(2));


%% check the effect of some distortion on these metrics
fprintf('Separation metrics without distortion\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30s(1), sir30s(1), sar30s(1));
%save orig separation for later
origsep = nmfrec;

%% distortion one: set 1 sample to 0 in each frame
nmfrec = origsep;
mask = ones(size(nmfrec{1}));
positions = [1:hp:totlength];
mask(positions) = 0;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) .* mask;
% soundsc(nmfrec{1}, 44100)
[sdr30d1, sir30d1, sar30d1, perm30d1] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with 1/hopsize sample removal\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d1(1), sir30d1(1), sar30d1(1));
% compute PEASS in one go
[ sdrP2d1, sirP2d1, sarP2d1, isrP2d1, ~, ~, ~, ~, OPSd1, TPSd1, IPSd1, APSd1 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);


%% now double the distortion
nmfrec = origsep;
mask = ones(size(nmfrec{1}));
positions = [1:hp/2:totlength];
mask(positions) = 0;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) .* mask;
[sdr30d12, sir30d12, sar30d12, perm30d12] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with 2/hopsize sample removal\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d12(1), sir30d12(1), sar30d12(1));
% compute PEASS in one go
[ sdrP2d12, sirP2d12, sarP2d12, isrP2d12, ~, ~, ~, ~, OPSd12, TPSd12, IPSd12, APSd12 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% now put the 0 samples next to eachother
nmfrec = origsep;
mask = ones(size(nmfrec{1}));
positions = [1:hp:totlength];
mask(positions) = 0; mask(positions+1) = 0;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) .* mask;
[sdr30d13, sir30d13, sar30d13, perm30d13] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with 2/hopsize consecutive sample removal\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d13(1), sir30d13(1), sar30d13(1));
% compute PEASS in one go
[ sdrP2d13, sirP2d13, sarP2d13, isrP2d13, ~, ~, ~, ~, OPSd13, TPSd13, IPSd13, APSd13 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion two: blend with white noise
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
nmfrec{1} = nmfrec{1} + rand(size(nmfrec{1})).*sqrt(mean(nmfrec{1}.^2)); % = RMS
%rescale for PEASS
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
[sdr30d2, sir30d2, sar30d2, perm30d2] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with 1xRMS white noise added\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d2(1), sir30d2(1), sar30d2(1));
% compute PEASS in one go
[ sdrP2d2, sirP2d2, sarP2d2, isrP2d2, ~, ~, ~, ~, OPSd2, TPSd2, IPSd2, APSd2 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% white noise at twice the RMS
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
nmfrec{1} = nmfrec{1} + rand(size(nmfrec{1})).*sqrt(mean(nmfrec{1}.^2)).*2; % = RMS
%rescale for PEASS
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
[sdr30d22, sir30d22, sar30d22, perm30d22] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with 2xRMS white noise added\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d22(1), sir30d22(1), sar30d22(1));
% compute PEASS in one go
[ sdrP2d22, sirP2d22, sarP2d22, isrP2d22, ~, ~, ~, ~, OPSd22, TPSd22, IPSd22, APSd22 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion three: add sinusoid at medium freq
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) + sin(2*pi*440*[1:1:totlength]./44100);
%rescale for PEASS
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
[sdr30d3, sir30d3, sar30d3, perm30d3] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with medium sinusoid tone added\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d3(1), sir30d3(1), sar30d3(1));
% compute PEASS in one go
[ sdrP2d3, sirP2d3, sarP2d3, isrP2d3, ~, ~, ~, ~, OPSd3, TPSd3, IPSd3, APSd3 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion four: add sinusoid at very-low freq
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) + sin(2*pi*20*[1:1:totlength]./44100);
%rescale for PEASS
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
[sdr30d4, sir30d4, sar30d4, perm30d4] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with very low sinusoid tone added\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d4(1), sir30d4(1), sar30d4(1));
% compute PEASS in one go
[ sdrP2d4, sirP2d4, sarP2d4, isrP2d4, ~, ~, ~, ~, OPSd4, TPSd4, IPSd4, APSd4 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion five: add sinusoid at ultra-low freq
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) + sin(2*pi*1*[1:1:totlength]./44100);
%rescale for PEASS
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1}));
[sdr30d5, sir30d5, sar30d5, perm30d5] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with ultra low sinusoid tone added\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d5(1), sir30d5(1), sar30d5(1));
% compute PEASS in one go
[ sdrP2d5, sirP2d5, sarP2d5, isrP2d5, ~, ~, ~, ~, OPSd5, TPSd5, IPSd5, APSd5 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion six: modulate with sinusoid at medium freq
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) .* sin(2*pi*440*[1:1:totlength]./44100);
[sdr30d6, sir30d6, sar30d6, perm30d6] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics modulated with medium sinusoid tone\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d6(1), sir30d6(1), sar30d6(1));
% compute PEASS in one go
[ sdrP2d6, sirP2d6, sarP2d6, isrP2d6, ~, ~, ~, ~, OPSd6, TPSd6, IPSd6, APSd6 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion seven: modulate sinusoid at very-low freq
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) .* sin(2*pi*20*[1:1:totlength]./44100);
[sdr30d7, sir30d7, sar30d7, perm30d7] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics modulated with very low sinusoid tone\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d7(1), sir30d7(1), sar30d7(1));
% compute PEASS in one go
[ sdrP2d7, sirP2d7, sarP2d7, isrP2d7, ~, ~, ~, ~, OPSd7, TPSd7, IPSd7, APSd7 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion eight: modulate sinusoid at ultra-low freq
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})) .* sin(2*pi*1*[1:1:totlength]./44100);
[sdr30d8, sir30d8, sar30d8, perm30d8] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics modulated with ultra low sinusoid tone\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d8(1), sir30d8(1), sar30d8(1));
% compute PEASS in one go
[ sdrP2d8, sirP2d8, sarP2d8, isrP2d8, ~, ~, ~, ~, OPSd8, TPSd8, IPSd8, APSd8 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion nine: add a DC offset
nmfrec = origsep;
% avoid clipping: scale to 0.9 then add 0.1
nmfrec{1} = nmfrec{1}./max(abs(nmfrec{1})).*0.9 + 0.1;
[sdr30d9, sir30d9, sar30d9, perm30d9] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics with +0.1 DC offset\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30d9(1), sir30d9(1), sar30d9(1));
% compute PEASS in one go
[ sdrP2d9, sirP2d9, sarP2d9, isrP2d9, ~, ~, ~, ~, OPSd9, TPSd9, IPSd9, APSd9 ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);

%% distortion ten: different gain
nmfrec = origsep;
nmfrec{1} = nmfrec{1}./10;
nmfrec{2} = nmfrec{2}./10;
[sdr30s, sir30s, sar30s, perm30s] =  bss_eval_sources(cell2mat(nmfrec'), cell2mat(arpegsrc'));
fprintf('Separation metrics at 1/10 scale.\n');
fprintf('BSS30s of source 1: SDR %f , SIR %f, SAR %f\n', sdr30s(1), sir30s(1), sar30s(1));
% compute PEASS in one go
[ sdrP2s, sirP2s, sarP2s, isrP2s, ~, ~, ~, ~, OPSs, TPSs, IPSs, APSs ] = EvalInfSepPEASS(nmfrec, arpegsrc, options);


% soundsc(nmfrec{1}, 44100)

%% print for BSS_EVAL LaTeX table
fprintf('1. Original separation result & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30s(1), sir30s(1), sar30s(1));
fprintf('2. Put 1 in 256 samples to 0 (crackle noise) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d1(1), sir30d1(1), sar30d1(1));
fprintf('3. Put 1 in 128 samples to 0 (more crackle noise) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d12(1), sir30d12(1), sar30d12(1));
fprintf('4. Add white noise scaled to 1 time the RMS of the signal & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d2(1), sir30d2(1), sar30d2(1));
fprintf('5. Add white noise scaled to twice the RMS of the signal & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d22(1), sir30d22(1), sar30d22(1));
fprintf('6. Add 440Hz sinusoid to signal (loud beep tone) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d3(1), sir30d3(1), sar30d3(1));
fprintf('7. Add 20Hz sinusoid to signal (barely audible rumble) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d4(1), sir30d4(1), sar30d4(1));
fprintf('8. Add 1Hz sinusoid to signal (not audible) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d5(1), sir30d5(1), sar30d5(1));
fprintf('9. Modulation with 440Hz sinusoid (bell-like resonance) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d6(1), sir30d6(1), sar30d6(1));
fprintf('10. Modulation with 20Hz sinusoid (``flutter'''' effect) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d7(1), sir30d7(1), sar30d7(1));
fprintf('11. Modulation with 1Hz sinusoid (loudness ``wobble\'''') & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d8(1), sir30d8(1), sar30d8(1));
fprintf('12. Add +0.1 DC offset to signal (not audible) & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30d9(1), sir30d9(1), sar30d9(1));
fprintf('13. Scaled to 1/10 of original gain & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdr30s(1), sir30s(1), sar30s(1));


%% print for PEASS latex table

fprintf('1. (Original) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2(1), sirP2(1), sarP2(1), isrP2(1), OPS(1), TPS(1), IPS(1), APS(1));
fprintf('2. (crackle, 1/256 sample) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d1(1), sirP2d1(1), sarP2d1(1), isrP2d1(1), OPSd1(1), TPSd1(1), IPSd1(1), APSd1(1));
fprintf('3. (crackle, 1/128 sample) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d12(1), sirP2d12(1), sarP2d12(1), isrP2d12(1), OPSd12(1), TPSd12(1), IPSd12(1), APSd12(1));
fprintf('4. (white noise, 1xRMS) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d2(1), sirP2d2(1), sarP2d2(1), isrP2d2(1), OPSd2(1), TPSd2(1), IPSd2(1), APSd2(1));
fprintf('5. (white noise, 2xRMS) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d22(1), sirP2d22(1), sarP2d22(1), isrP2d22(1), OPSd22(1), TPSd22(1), IPSd22(1), APSd22(1));
fprintf('6. (Add 440Hz tone) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d3(1), sirP2d3(1), sarP2d3(1), isrP2d3(1), OPSd3(1), TPSd3(1), IPSd3(1), APSd3(1));
fprintf('7. (Add 20Hz tone) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d4(1), sirP2d4(1), sarP2d4(1), isrP2d4(1), OPSd4(1), TPSd4(1), IPSd4(1), APSd4(1));
fprintf('8. (Add 1Hz tone) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d5(1), sirP2d5(1), sarP2d5(1), isrP2d5(1), OPSd5(1), TPSd5(1), IPSd5(1), APSd5(1));
fprintf('9. (Modulate with 440Hz) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d6(1), sirP2d6(1), sarP2d6(1), isrP2d6(1), OPSd6(1), TPSd6(1), IPSd6(1), APSd6(1));
fprintf('10. (Modulate with 20Hz) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d7(1), sirP2d7(1), sarP2d7(1), isrP2d7(1), OPSd7(1), TPSd7(1), IPSd7(1), APSd7(1));
fprintf('11. (Modulate with 1Hz) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d8(1), sirP2d8(1), sarP2d8(1), isrP2d8(1), OPSd8(1), TPSd8(1), IPSd8(1), APSd8(1));
fprintf('12. (Add +0.1 DC offset) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2d9(1), sirP2d9(1), sarP2d9(1), isrP2d9(1), OPSd9(1), TPSd9(1), IPSd9(1), APSd9(1));
fprintf('13. (Scaled to 1/10th) & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ & $%2.2f$ \\\\\n\\hline\n', sdrP2s(1), sirP2s(1), sarP2s(1), isrP2s(1), OPSs(1), TPSs(1), IPSs(1), APSs(1));


