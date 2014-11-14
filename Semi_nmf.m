function [sum_err,num_ite,B,C] = Semi_nmf(A,k,itor,tor_err)
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
for p=1:k
    B(:,p)=B(:,p)/sqrt(Num_k(p));
    %sqrt(Num_k(p))
end
B;
B = B+0.2;
C = A*B/(B'*B);
Numerator_Positive = (abs(A'*C)+A'*C)/2;
Numerator_Negative = B*(abs(C'*C)-C'*C)/2;
Denominator_Positive = B*(abs(C'*C)+C'*C)/2;
Denominator_Negative = (abs(A'*C)-A'*C)/2;
B = B.*sqrt((Numerator_Negative+Numerator_Positive)./(Denominator_Negative+Denominator_Positive));
Err = zeros(1,itor);
Err(1) = sum(sum((A-C*B').^2));
for j=2:itor
    
    C = A*B/(B'*B);
    Numerator_Positive = (abs(A'*C)+A'*C)/2;
    Numerator_Negative = B*(abs(C'*C)-C'*C)/2;
    Denominator_Positive = B*(abs(C'*C)+C'*C)/2;
    Denominator_Negative = (abs(A'*C)-A'*C)/2;
    B = B.*sqrt((Numerator_Negative+Numerator_Positive)./(Denominator_Negative+Denominator_Positive));

   %B = B.*((A'*C)./(B*B'*A'*C));
   %C = C.*((A*B)./(C*B'*B));
   Err(j) = sum(sum((A-C*B').^2));
   sum_err=Err(j);
   num_ite=j;
   if (Err(j-1)-Err(j))<tor_err
       break;
   end
end