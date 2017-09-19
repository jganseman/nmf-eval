function wavwrite(y,Fs,filename)
% Function wavwrite provided for forward compatibility with Matlab versions
% after 2015 
% (wavwrite was removed in 2015b)
% (audiowrite was introduced in 2012b)

% author: Joachim Ganseman

audiowrite(filename, y, Fs);