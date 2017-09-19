function [output, sampleRate] = wavread(file, varargin)
% Function wavread provided for forward compatibility with Matlab versions
% after 2015 
% (wavread was removed in 2015b)
% (audioread was introduced in 2012b)

% author: Joachim Ganseman

% special case: if varargin is 'size', return [samples channels] of audiodata
if ((length(varargin) == 1) && strcmp(varargin{1},'size'))
    info = audioinfo(file);
    output = [ info.TotalSamples/info.NumChannels, info.NumChannels] ; 
    %TODO necessary to divide by NumChannels?
else
    [output, sampleRate] = audioread(file, varargin{:});
end

