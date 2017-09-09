function [xr,win_pos] =  istft (X,hopfac,winlen,type)
% xr =  istft (X,hopfac,winlen)
% this function calculate the inverse STFT for a STFT matrix
% X - STFT matrix (bins 1:nfft/2+1)
% hopfac - hop factor. This is an integer specifying the number of analysis hops
% occurring within a signal frame that is one window in length. In other words,
% winlen/hopfac is the hop size in saamples and winlen*(1-1/hopfac) is the overlap
% in samples.
% winlen - window length in samples
% type - 'perfect' or 'smooth'
% (c) Shlomo Dubnov sdubnov@ucsd.edu

% Addition by Joachim Ganseman (Dec. 2013): 
% - This code has been corrected to do a unitary transform in the 'smooth'
% case as well. See explanataion at the bottom of this file.
% - Updates: 'hanning' -> 'hann' , and corrected string comparisons

X = [X; conj(X(end-1:-1:2,:))];

if nargin < 2,
    hopfac = 2;
end
if nargin < 3,
    winlen = size(X,1);
end
if nargin < 4,
    type = 'perfect';
end

hop = winlen/hopfac;
bmat = ifft(X);
%STFT = real(ifft(X));

[M N] = size(bmat);
nfft = M;

xr = zeros (1,N*hop + nfft);
win_pos = [1: hop: length(xr) - nfft];
if strcmp(type, 'perfect'),
    win = ones(winlen,1);
elseif strcmp(type, 'smooth');
    win = hann(winlen,'periodic'); %second smoothing window
else
    error('no such istft type')
end


for i=1:length(win_pos)
    xr(win_pos(i):win_pos(i)+nfft-1) = xr(win_pos(i):win_pos(i)+nfft-1) + bmat(:,i)'.*win';
end

% NOTE by Joachim Ganseman:
% " xr = real(xr)/hopfac*2; " works in case 'perfect' is chosen, because
% the '2' corrects for the application of the window in the forward transform, (by
% multiplying with the average value of the hann window (0.5). However, it is wrong
% in case of 'smooth' transform because then the unwindowing step at the
% end is not taken into account.

% Therefore, generalize this: apply the hop factor correction,
% and another normalization based on the window correction factors in both
% forward and backward transforms to obtain unitarity

if strcmp(type, 'perfect'),
    xr = real(xr)/hopfac*2;
elseif strcmp(type, 'smooth');
    xr = real(xr) / (hopfac*mean(win.^2));
end

% TODO: this should be generalized, so as to work with windows in general
% However, when different windows are used at analysis and synthesis, this
% isn't all too easy.


