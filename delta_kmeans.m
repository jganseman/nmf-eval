function [sum_err,num_ite,B,C] = delta_kmeans(A,k,itor,tor_err)
[IDX,C]=kmeans(A',k);
B = zeros(size(A,1),k);
for i=1:size(IDX)
    B(i,IDX(i))=1;
end
B = B';
C = C';
Err = zeros(1,itor);
Err(1) = sum(sum((A-C*B).^2));
for j=2:itor
   B=B.*((C'*A)./(C'*C*B));
   C=C.*((A*B')./(C*B*B'));
   Err(j) = sum(sum((A-C*B).^2));
   sum_err=Err(j);
   num_ite=j;
   if (Err(j-1)-Err(j))<tor_err
       break;
   end
end
    