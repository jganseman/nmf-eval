
% add paths
addpath ./scatt_functions
addpath ./bss_eval

% here goes a link to the SPAMS library
addpath ../spams-matlab/build/

% define scattering parameters and construct filters
eps=1e-3;
fs = 16000;
Npad = 2^15;
T = 2048;
scparam.N = Npad;
scparam.T = T;
scparam.Q = 32;
filts = create_scattfilters( scparam );


% Load precomputed dictionaries go to precompute_dicts.m to see how this is done (SPAMS needed)
% These files contain the the dictionaries for both scattering levels and the parametrs used for computing them (param1 and param2)
load ('pretrained_dics/dictionaries_NMF_scatt2_s7')
Dnmf11 = Dnmf1;
Dnmf21 = Dnmf2;
clear Dnmf1 Dnmf2


load ('pretrained_dics/dictionaries_NMF_scatt2_s1')
Dnmf12 = Dnmf1;
Dnmf22 = Dnmf2;
clear Dnmf1 Dnmf2

% load renormalization parameters (see precompute_dicts.m )
load('pretrained_dics/renorm_params')

% use this pair of files for an example of good performance
file_1 = 'audio_files/sbba2n.wav'
file_2 = 'audio_files/sbbh2n.wav'

% use this pair of files for a more challenging case
%file_1 = 'audio_files/bgat5a.wav'
%file_2 = 'audio_files/bgbh5s.wav'

% Load files and resample to fs = 16kHz
[x1, Fs] = audioread( file_1  );
x1 = resample(x1,fs,Fs);
x1 = x1(:)'; T1 = length(x1 );

[x2, Fs] = audioread( file_2 );
x2 = resample(x2,fs,Fs);
x2 = x2(:)'; T2 = length(x2);


% create the mixture
T = min([T1,T2,Npad]);
x1 = x1(1:T);
x2 = x2(1:T);
mix = (x1+x2);

% perform separation
[speech1, speech2, xest1, xest2] = demix_scatt2top(mix, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, eps, filts, scparam, param1, param2, Npad);


% evaluate results
Parms_scatt2  =  BSS_EVAL(x1', x2', speech1(1:T)', speech2(1:T)', mix');
Parms_scatt1 =  BSS_EVAL(x1', x2', xest1(1:T)', xest2(1:T)', mix');

Parms_scatt2
Parms_scatt1



    
    

