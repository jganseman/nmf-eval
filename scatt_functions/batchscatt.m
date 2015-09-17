function [out,pout] = batchscatt(in, filts, options)

%first order scattering with oversampling

J=size(filts.psi,2);
os = getoptions(options,'os',1);
Q = getoptions(options,'Q',32);
T = getoptions(options,'T',4096);
dse = getoptions(options,'dse',round(T*2^(-floor(J/Q)-os)));

[N, L]=size(in);

Neff = floor(N/dse);
Xf = fft(in);
out=zeros(J+1,Neff, L);
pout=zeros(J+1,N, L);
for j=1:J
tmp = (ifft(Xf.*repmat(filts.psi{j},1,L)));
pout(j,:,:)=exp(i*angle(tmp));
tmpds = abs(tmp(1:dse:end,:));
out(j,:,:)=tmpds;
end
tmp = (ifft(Xf.*repmat(filts.phi,1,L)));
pout(J+1,:,:)=exp(i*angle(tmp));
tmpds = abs(tmp(1:dse:end,:));
out(J+1,:,:)=tmpds;




%Neff = floor(N/T);
%
%Xf = fft(in);
%ncoeffs=0;
%for j=1:J
%idse = 2^(floor((J-j)/Q));
%ncoeffs = ncoeffs + idse;
%end
%out=zeros(ncoeffs,Neff, L);
%
%rast=1;
%for j=1:J
%tmp = abs(ifft(Xf.*repmat(filts.psi{j},1,L)));
%idse = 2^(floor((J-j)/Q));
%dse = T/idse; %*2^((j-J)/Q));  
%tmpds = tmp(1:dse:end,:);
%p = size(tmpds,1)/Neff;
%out(rast:rast+p-1,:,:) = reshape(tmpds, p, Neff, L);
%rast = rast + p;
%end
%
%


