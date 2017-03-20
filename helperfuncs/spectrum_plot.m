function spectrum_plot(X,raise_factor)

% Plot the given magnitude spectrum using nice looking options 
%(axis xy and colormap(hot).
%
% spectrum_plot(X,raise_factor)
% 
% X = input magnitude spectrogram
% raise_factor = the power to which we raise the spectrogram (for better
% clarity) (default = 1)

if nargin < 2
    raise_factor = 1;
end

colormap(1-gray)

imagesc(X.^raise_factor)
%edit by JOA: use smaragdis' scaling method from cplca.m:
%imagesc( reshape( X, size( X, 1), [])./max(X(:)))
%

axis xy
% colormap(hot)         % adapted by JOA to:
%colormap(1-contrast(X))         % use this when plotting log spectrograms
% also put colormap first, because imagesc is based on current colormap