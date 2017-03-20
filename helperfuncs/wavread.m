function [output sampleRate] = wavread(file, varargin)
% Function wavread provided for forward compatibility with Matlab versions
% after 2015 
% (wavread was removed in 2015b)
% (audioread was introduced in 2012b)

% author: Joachim Ganseman

[output sampleRate] = audioread(file, varargin{:});
