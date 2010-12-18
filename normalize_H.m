function H = normalize_H(H, type)
% function H = normalize_H(H, type)
%
% Normalize H using type which can be:
%  -1   - use 1-norm of vectorized matrix
%  -2   - use 2-norm of vectorized matrix
%   1   - use 1-norm of each row [default]
%   2   - use 2-norm or each row
%   k   - multiply the 1-norm of each row by k
%   'a' - means make sum(H(:))=1 (same as -1)
%
% 2010-01-14 Graham Grindlay (grindlay@ee.columbia.edu)

% Copyright (C) 2008-2028 Graham Grindlay (grindlay@ee.columbia.edu)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin < 2
    type = 0;
end

switch type
    case -1
        H = H ./ norm(H(:),1);
        
    case -2
        H = H ./ norm(H(:),2);
        
    case 1
        for i = 1:size(H,1)
            H(i,:) = H(i,:) ./ norm(H(i,:),1);
        end
        
    case 2
        for i = 1:size(H,1)
            H(i,:) = H(i,:) ./ norm(H(i,:),2);
        end
        
    case 'a' % TODO: clear this all up
        H = H./sum(H(:));
        
    case 0
        
    otherwise % TODO: get rid of this
        for i = 1:size(H,1)
            H(i,:) = type*H(i,:) ./ norm(H(i,:),1);
        end
end
        