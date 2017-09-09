function [ output_args ] = getJGDB128( input_args )
%GETJGDB128 Returns the place of the data directory for JGDB 128 database
%   Use this function to define a path to the example data files

% Note: if necessary, replace the paths in this file with the paths
% that work on your local machine!

if ismac
    str = '/Volumes/JG-DATA/databases/jgdb128/';
elseif isunix
    str = '/data/databases/jgdb128/';
elseif ispc
    str = 'F:\databases\jgdb128\';
end

end

