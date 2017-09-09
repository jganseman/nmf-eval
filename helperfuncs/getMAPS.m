function str = getMAPS()
%GETMAPS Returns the place of the data directory for MAPS database
%   Use this function to define a path to the example data files

% Note: if necessary, replace the paths in this file with the paths
% that work on your local machine!

if ismac
    str = '/Volumes/JG-DATA/databases/MAPS/';
elseif isunix
    str = '/data/databases/MAPS/';
elseif ispc
    str = 'F:\';
end

end

