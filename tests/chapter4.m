% Chapter 4: creating the images for Chapter 4 of my PhD Thesis.
% (c) Joachim Ganseman, september 2017

% NOTE: AFTER RUNNING, EXTERNALIZE THE DATA TABLE READING in 2 steps:
% 1. add the following code before subimporting the created tex file:
%     \pgfplotstableread{./images/chapter2/audiosignalexcerpt-1.tsv}\loadedtable
% 2. in the tex file itself, change the following:
%     table[row sep=crcr,format=inline]{\loadedtable}; 

clear;

%%%% PART ONE: illustrate audio-to-score alignment %%%%
% This is based on [Ganseman 2010 ICMC]

PRINTTOFILE=1;

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
filename='../../../thesis/images/chapter4/dtwexample.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

