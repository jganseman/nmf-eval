% Short-Time Fourier Transform (STFT)
%     [X_,m_] = STFT(x_,w_,R,M) produces the STFT X_ of vector x_, where w_
%     is the analysis window, R the hop size, and M, the FFT length. m_ is
%     the vector of time indices corresponding to the first sample of each
%     short-time spectrum.

%     The function was written to support multi-channel input, the number
%     of columns in x_ being interpreted as the number of channels.
%     The size of X_ will therefore be (M,chans,U), where chans is the
%     number of channels, and U, the minimal number of windows for the
%     transparent reconstruction of x_ upon inverse STFT - provided that
%     w_ and R were chosen in a way to make transparent reconstruction
%     possible.

% Adaptation by Joachim Ganseman (Dec. 2013): cut redundant part of
% spectrogam at the end. Also, comments, and removed the unnessecary M/2pi

function [X_,m_] = STFT(x_,w_,R,M)

[in_N,chans] = size(x_);
ww_ = repmat(w_,1,chans);
N = length(w_);
if nargin<4, M = 2^ceil(log2(N)); end

O = N/R;                    % hops per frame
if mod(O,1), error('Window length must be multiple of hop size'); end
U = ceil(in_N/R)+O-1;       % necessary nr of frames
u_ = (1:U)';
m_ = (u_-O)*R;              %vector of indices, [-768 -512 ... siglength]
l_Z = (O-1)*R;              % nr samples to left pad
r_Z = (O+floor(in_N/R))*R-in_N;     % nr samples to right pad

zxz_ = [zeros(l_Z,chans);x_;zeros(r_Z,chans)];  % padded signal
ZNZ = in_N+l_Z+r_Z;         % size of padded signal

%M2pi = M/2/pi;
X_ = zeros(M,chans,U);
for u=u_'
    n_ = m_(u)+(0:N-1)';            % vector of M sample nrs
    N_x_ = zxz_(n_+l_Z+1,:);        % cut them from padded signal
    X_(:,:,u) = fft(N_x_.*ww_,M);   % perform fft
end

% JGA addition: cut redundant part of spectrogram
X_ = X_(1:M/2+1,:,:);