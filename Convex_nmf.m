function [sum_err,num_ite,W,G,F] = Convex_nmf(A,k,itor,tor_err)
[IDX,C]=kmeans(A',k);
Num_k=zeros(1,k);
for i = 1:k
    for j = 1:size(IDX)
        if i == IDX(j)
          Num_k(i)=Num_k(i)+1;
        end
    end
end

% for i = 1:k
%     Num_k(i)
% end

B = zeros(size(A,2),k);

for i=1:size(IDX)
    B(i,IDX(i))=1;
end
B;

G=B+0.2;
D=zeros(k,k);
for i=1:k
    D(i,i)=Num_k(i);
end

W = G/D;
Err = zeros(1,itor);
Err(1) = sum(sum((A-A*W*G').^2));
for j=2:itor
   
    Positive = (abs(A'*A)+A'*A)/2;
    
    Negative = (abs(A'*A)-A'*A)/2;
    
    G = G.*sqrt((Positive*W+G*W'*Negative*W)./(Negative*W+G*W'*Positive*W));
    W = W.*sqrt((Positive*G+Negative*W*G'*G)./(Negative*G+Positive*W*G'*G));
    F = A*W;
   Err(j) = sum(sum((A-A*W*G').^2));
   sum_err=Err(j);
   num_ite=j;
   if (Err(j-1)-Err(j))<tor_err
       break;
   end
end