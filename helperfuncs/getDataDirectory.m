function str = getDataDirectory()
%GETDATADIRECTORY Returns the place of the data directory
%   Use this function to define a path to the example data files

if ismac
    str = './testdata/';
elseif isunix
    str = './testdata/';
elseif ispc
    str = '.\testdata\';
end

% if at C4DM, use str = '/Volumes/c4dm-scratch/jga/data/';
% if at C4DM on linux, use str = '/c4dm-scratch/jga/data/';
% if at home, use str = '/Users/jg/Dropbox/data/'; 
% if at home on pc, use, str = 'C:\Users\JG\Dropbox\data\';
    
%end

end

