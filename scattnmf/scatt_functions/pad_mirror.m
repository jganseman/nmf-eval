function out=pad_mirror(in, L)


out=zeros(L,1);
T=length(in);
if T > L
out = in(1+floor((T-L)/2):floor((T-L)/2)+L);
else

boo=linspace(in(end),in(1),L-T);
out(1:T)=in;
out(T+1:end)=boo;

end


