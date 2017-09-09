function [ output_args ] = getNSynth( input_args )
%GETNSYNTH Returns the place of the data directory for NSynth database
%   Use this function to define a path to the example data files

% Note: if necessary, replace the paths in this file with the paths
% that work on your local machine!

if ismac
    str = '/Volumes/JG-DATA/databases/nsynth/';
elseif isunix
    str = '/data/databases/nsynth/';
elseif ispc
    str = 'F:\databases\nsynth\';
end

end

