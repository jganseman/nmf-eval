function [sum_err,num_ite,F,S,G] = Bi_orthogonal(A,k,l,itor,tor_err)
%
%
%
%Implement Bi_Orthogonal tri-factorization, see Chris Ding et al Orthgonal
%Nonnegative Matrix Tri-factorizations for Clustering,2006.
%Object  Min||A-FSG'|| s.t. F,S,G>=0, G'*G=I,F'*F=I
%
%
%
%Inputs:
%   A     [mat]  - Input matrix (n by m)   
%   k     [num]  - Rank of Decomposition of column clustering of A
%   l     [num]  - Rank of Decomposition of Row clustering of A
%   itor  [num]  - Max times of iteration for the update 
%   tor_err  [num]  - Number to determine convergence, if the difference of
%                     two consective lost function is less than this, stop
%                     the update
%
%
%Outputs:
%   sum_err     [num]  - value of the lost function||A-FSG'||
%   num_ite     [num]  - times of update
%   F     [mat]  - Output matrix (n by l) indocator membership for the row
%                  clustering of A
%   G     [mat]  - Output matrix (m by k) indocator membership for the col
%                  -umn clustering of A
%   S     [mat]  - Output matrix (l by k) rescale factor
%
%
%
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


% do the kmeans clustering for the column of A
%Num_Column_k denotes the numbers for each clustering
IDX_Column=kmeans(A',k);
Num_Column_k=zeros(1,k);
for i = 1:k
    for j = 1:size(IDX_Column)
        if i == IDX_Column(j)
          Num_Column_k(i)=Num_Column_k(i)+1;
        end
    end
end

% for i = 1:k
%     Num_k(i)
% end


% transfer from the output IDX to G,set the accordingly element to 1 
% and to guarantee the orthogonality of G, divide each column by sqrt(nk)
% then add 0.2 to each element
G = zeros(size(A,2),k);

for i=1:size(IDX_Column)
    G(i,IDX_Column(i))=1;
end
G;
for p=1:k
    G(:,p)=G(:,p)/sqrt(Num_Column_k(p));
    %sqrt(Num_k(p))
end
G;
G = G+0.2;

%--------------------------------------------------------------------------
% the following code do the similar thing above, to get the row indicator 
% membership for A by doing the row kmeans clustering
%--------------------------------------------------------------------------
IDX_Row=kmeans(A,l);
Num_Row_k=zeros(1,l);
for i = 1:l
    for j = 1:size(IDX_Row)
        if i == IDX_Row(j)
          Num_Row_k(i)=Num_Row_k(i)+1;
        end
    end
end

% for i = 1:k
%     Num_k(i)
% end

F = zeros(size(A,1),l);

for i=1:size(IDX_Row)
    F(i,IDX_Row(i))=1;
end
F;
for p=1:l
    F(:,p)=F(:,p)/sqrt(Num_Row_k(p));
    %sqrt(Num_k(p))
end
F;
F = F+0.2;
%-----------------------------------------------------------
%initialize S
%-----------------------------------------------------------
S = F'*A*G;
%--------------------------------------------------------------------------
%Err donotes the lost function ||A-FSG'|| for each step of update
%do the update according to the algorithm provided in the paper, and when
%the difference of two consective lost function is less than tor_err, stop
%the update
%--------------------------------------------------------------------------

Err = zeros(1,itor);
Err(1) = sum(sum((A-F*S*G').^2));
for j=2:itor
   G = G.*((A'*F*S)./(G*G'*A'*F*S));
   F = F.*((A*G*S')./(F*F'*A*G*S'));
   S = S.*((F'*A*G)./(F'*F*S*G'*G));
   Err(j) = sum(sum((A-F*S*G').^2));
   sum_err=Err(j);
   num_ite=j;
   if (Err(j-1)-Err(j))<tor_err
       break;
   end
end
    