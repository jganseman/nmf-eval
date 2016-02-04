function out=haarinv(in)

Jhaar=size(in,2)-1;

N=size(in{1},1)+1;
L=size(in{1},2);

out=zeros(N,L);
out(2^(Jhaar-1):end-2^(Jhaar-1),:)=in{end};

for j=Jhaar:-1:1
nout=zeros(N,L);
tmp=0*nout;
tmp(2^(j-1):end-2^(j-1),:)=in{j};
I=1:ceil(2^(j-0));
nout(I,:)=out(I+floor(2^(j-2)),:)-tmp(I+floor(2^(j-2)),:);
I=ceil(2^(j-0))+1:N-floor(2^(j-0));
nout(I,:)=.5*(out(I+floor(2^(j-2)),:)-tmp(I+floor(2^(j-2)),:)+out(I-ceil(2^(j-2)),:)+tmp(I-ceil(2^(j-2)),:));
I=N-floor(2^(j-0)):N;
nout(I,:)=out(I-ceil(2^(j-2)),:)+tmp(I-ceil(2^(j-2)),:);
out = nout;
end



