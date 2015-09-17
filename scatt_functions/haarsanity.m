%%%% simple test to verify haar TI inversion
clear all;

options.null=0;

N=1024;

Jhaar=getoptions(options,'Jhaar',5);
for j=1:Jhaar
slice=zeros(2^j,1);
slice(1:2^(j-1))=2^(-j);
slice(2^(j-1)+1:end)=-2^(-j);
filts.haar{j}=slice;
filts.ds(j)=2^(j-1);
end
filts.haar{Jhaar+1}=2^(-Jhaar)*ones(2^Jhaar,1);
filts.ds(Jhaar+1)=2^(Jhaar-1);



in=randn(N,100);
%in=0*in;
%in(2)=1;
%in(end-1)=1;

%forward transform
for j=1:Jhaar+1
fin{j} = conv2(in,filts.haar{j},'valid');
end

pout=haarinv(fin);

if 0
%inverse transform
out=zeros(N,1);
%out(2^(Jhaar-1)+0:end-2^(Jhaar-1))=fin{end};

out(2^(Jhaar-1):end-2^(Jhaar-1))=fin{end};

for j=Jhaar:-1:1
nout=zeros(N,1);
tmp=0*nout;
tmp(2^(j-1):end-2^(j-1))=fin{j};
%I=2^(j-1)+0:2^(j)+2^(j-1);
I=1:ceil(2^(j-0));
%nout(I)=out(I+floor(2^(j-2)))+fin{j}(I+floor(2^(j-2))-2^(j-1) );
nout(I)=out(I+floor(2^(j-2)))-tmp(I+floor(2^(j-2)));
%I=2^(j)+2^(j-1)+1:N-2^(j)-2^(j-1);
I=ceil(2^(j-0))+1:N-floor(2^(j-0));
%nout(I)=.5*(out(I+floor(2^(j-2)))+fin{j}(I+floor(2^(j-2))-2^(j-1) )+out(I-ceil(2^(j-2)))-fin{j}(I-ceil(2^(j-2))-2^(j-1) ));
nout(I)=.5*(out(I+floor(2^(j-2))+0)-tmp(I+floor(2^(j-2))+0)+out(I-ceil(2^(j-2)))+tmp(I-ceil(2^(j-2)) ));
%I=N-2^(j-1)-2^(j)+1:N-2^(j-1)+1;
I=N-floor(2^(j-0)):N;%2^(j-1)-2^(j)+1:N-2^(j-1)+1;
%nout(I)=out(I-ceil(2^(j-2))+1)-fin{j}(I-ceil(2^(j-2)+1) );
nout(I)=out(I-ceil(2^(j-2)))+tmp(I-ceil(2^(j-2)));
out = nout;
end


end



