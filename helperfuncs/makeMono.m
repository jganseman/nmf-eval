function outwav = makeMono(inwav)
%MAKEMONO  take a sound and make it mono
%   outwav = makeMono(inwav) takes a sound that has just been read in with
%   wavread (thus, in column format), and makes it mono (returning a single
%   column). To avoid clipping, all channels are added then divided by the
%   nr of channels. The function does nothing when input is already mono.
%
% author: Joachim Ganseman
% precondition: inwav comes straight from wavread (column-wise data)
% postcondition: outwav is a single column vector

if size(inwav, 2) > 1
   outwav = sum(inwav, 2) ./ size(inwav, 2);
else
    outwav = inwav;
end

%TODO make it work on cells containing sound too?

