function [A,X,cost] = SNMF_mult(varargin,Dataoptions)
%  multiplicative rule for KL-NMF: Y=AX
% usage SNMF_mult({Y,A,X,max_it},{'threshold',t,'updateBasis',1,'costfn','kl','min_it',100})
% Nasser Mohammadiha January 2010

names =  {'threshold','updateBasis','costfn','min_it'};
options = {0.01,1,'kl',120};


Y=varargin{1};A=varargin{2};X=varargin{3};max_it=varargin{4};
for i=1:2:length(Dataoptions),
    idx = strmatch(Dataoptions{i},names,'exact');
    if isempty(idx), error(['Unknown Option : ' Dataoptions{i}]); end;
    options{idx} =  Dataoptions{i+1};
end;
impr_threshold=options{1};updateBasis=options{2};
costfn=options{3};
min_it=options{4};
imp=1;cost=[];
[F,T]=size(Y);


for l=1:max_it
    X=(X.*(A'*(Y./(A*X+eps)))./(max(repmat(A'*ones(F,1),1,T),1e-3)));
    
    if(updateBasis == 1)
        A = (A.*((Y./(A*X + eps))*X')./(max(repmat(sum(X,2)',F,1),1e-3)));
        A = A*diag(1./(sum(A,1) + eps)); %normalize A such that each basis vector has a norm1=1
    end
    
    er_matrix=Y.*log((Y+eps)./(A*X+eps))-Y+A*X;
    cc(l)=sum(er_matrix(:));
    cost(l)=cc(l);
    if(imp<impr_threshold)
        if(l>min_it)
            break;
        end
    end
    if l>=2
        imp=abs(cost(l)-cost(l-1))/abs(cost(l-1));
    end
end

end