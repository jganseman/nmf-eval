function f = stft( x, sz, hp, pd, w, ll)
% Short-time Fourier transform
%
% function f = stft( x, sz, hp, pd, w, ll)
%
% Inputs: 
%  x   input time series (must be row vector), or input complex spectrogram (DC to Nyquist)
%  sz  size of the FFT
%  hp  hop size in samples
%  pd  pad size in samples
%  w   window to use (function name of data vector)
%  ll  highest frequency index to return
%
% Output:
%  f   complex STFT output (only DC to Nyquist components), or time series resynthesis

% Paris Smaragdis 2006-2008, paris@media.mit.edu

% Bugfix By Joachim Ganseman (Dec. 2013): windows should be periodic, as to
% ascertain the Constant-Overlap-Add (COLA) condition for perfect STFT
% reconstruction. This is an additional flag in the window functions.

% Addition by Joachim Ganseman (Dec. 2013): make the transform unitary.
% scale the output by the appropriate window correction factor, which is
% the RMS of the window: sqrt(mean(window^2)).

% Forward transform
if isreal( x)

	% Defaults
	if nargin < 5
		w = 1;
	end
	if nargin < 4
		pd = 0;
	end
	if nargin < 3
		hp = sz/2;
	end

	% Specified window is a string
	if isstr( w)
		w = feval( w, sz , 'periodic');      %JGA: periodic added.
	end

	% Orient and zero pad input
	if size( x, 1) > 1
		x = x';
	end
	x = [x zeros( 1, ceil( length(x)/sz)*sz-length(x))];
%	x = [zeros( 1, sz+pd) x zeros( 1, sz+pd)];

	% Pack frames into matrix
	if isa( x, 'single')
		s = zeros( sz, (length(x)-sz)/hp, 'single');
	else
		s = zeros( sz, (length(x)-sz)/hp);
	end
	j = 1;
	for i = sz:hp:length( x)
		s(:,j) = w .* x((i-sz+1):i).';
		j = j + 1;
	end

	% FFT it
	f = fft( s, sz+pd);

	% Chop redundant part
	f = f(1:end/2+1,:);

	% Chop again to given limits
	if nargin == 6
		f = f(1:ll,:);		
    end
    
    %JGA: add scaling by window correction factor
    f = f ./ sqrt(mean(w.^2));
	
	% Just plot
	if nargout == 0
		imagesc( log( abs(f))), axis xy
%		xlabel( 'Time (sec)')
%		ylabel( 'Frequency (Rad)')
%		set( gca, 'xtick
	end

% Inverse transform
else

	% Defaults
	if nargin < 5
		w = 1;
	end
	if nargin < 4
		pd = 0;
	end
	if nargin < 3
		hp = sz/2;
	end

	% Specified window is a string
	if isstr( w)
		w = feval( w, sz ,'periodic');      %JGA: periodic added.
	end

	% Ignore padded part
	if length( w) == sz
		w = [w; zeros( pd, 1)];
	end

	% Overlap add/window/replace conjugate part
	if isa( x, 'single')
		f = zeros( 1, (size(x,2)-1)*hp+sz+pd, 'single');
	else
		f = zeros( 1, (size(x,2)-1)*hp+sz+pd);
	end
	v = 1:sz+pd;
	for i = 1:size( x,2)
		f((i-1)*hp+v) = f((i-1)*hp+v) + ...
			(w .* real( ifft( [x(:,i); conj( x(end-1:-1:2,i))])))';
	end

	% Norm for overlap
	f = f / (sz/hp);
%	f = f(sz+pd+1:end-sz-2*pd);

    %JGA: add scaling by window correction factor
    f = f ./ sqrt(mean(w.^2));

end
