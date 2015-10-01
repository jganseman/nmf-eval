function [U,P] = fwdcqt(in, options, filts)

J1 = size(filts.psi,2);
[N, L]=size(in);
Xf = fft(in);

P=zeros(J1+1,N, L);

%first layer
for j=1:J1
tmp = (ifft(Xf.*repmat(filts.psi{j},1,L)));
scratch.X{1}(j,:,:) = tmp(1:filts.ds(j):end,:);
P(j,:,:)=exp(i*angle(tmp));
end
tmp = ifft(Xf.*repmat(filts.phi,1,L));
scratch.X{1}(J1+1,:,:) = tmp(1:filts.ds(J1+1):end,:);
P(J1+1,:,:)=exp(i*angle(tmp));

scratch.P = P;
U=abs(scratch.X{1});
U(end,:,:) = scratch.X{1}(end,:,:);

