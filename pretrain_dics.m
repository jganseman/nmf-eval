


    
representation = '/misc/vlgscratch3/LecunGroup/bruna/grid_data/scatt2_fs16_NFFT2048/';

%compute renormalization parameters    
if ~exist('stds1','var')
stds1 = std(X1,0,2);
stds2 = std(X2,0,2);
end


% this code assumes that data contains in the fields X1 and X2 the first and second level scattering coefficients for some training data.  
% these coefficients can be computed as:
% [X2, X1] = audioscatt_fwd_haar(pad_mirror(x',Npad), filts, options);
% where x is the audio signal in the time domain

Npad = 2^15;
    
options.renorm=1;
if options.renorm
   %renormalize data: whiten each frequency component.
   eps=2e-3;
   data.X1 = renorm_spect_data(data.X1, stds1, eps);
        
   eps=1e-3;
   data.X2 = renorm_spect_data(data.X2, stds2, eps);
end
   
    
%% train models
model = 'NMF-scatt2';
    
    
%%%%Plain NMF%%%%%%%
KK1 = [160];
LL1 = [0.1];
param1.K = KK1;
param1.posAlpha = 1;
param1.posD = 1;
param1.pos = 1;
param1.lambda = LL1;
param1.iter = 4000;
param1.numThreads=16;
param1.batchsize=512;

% need spams to compute this dictionary, but any algorithm for NMF would do     
Dnmf1 = mexTrainDL(abs(data.X1),param1);
    
KK2 = [768];
LL2 = [0.1];
param2.K = KK2;
param2.posAlpha = 1;
param2.posD = 1;
param2.pos = 1;
param2.lambda = LL2;
param2.iter = 4000;
param2.numThreads=16;
param2.batchsize=512;
    
% need spams to compute this dictionary, but any algorithm for NMF would do     
Dnmf2 = mexTrainDL(abs(data.X2),param2);
    


