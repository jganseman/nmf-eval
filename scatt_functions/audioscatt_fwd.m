function [S,scratch] = audioscatt_fwd(in, filts, options)
%this function computes a two-layer speech scattering transform 

%filts{1}: first layer CQT filters
%filts{2}: second layer across-octave filters
%filts{3}: second layer within-octave filters
%filts{4}: second layer dyadic temporal filters

% T : final downsampling factor
% dse : first order downsampling 
J1 = size(filts{1}.psi,2);
%os = getoptions(options,'os',1);
%T = getoptions(options,'T',4096);
%dse = getoptions(options,'dse',round(T*2^(-floor(J/Q)-os)));

[N, L]=size(in);
%Neff = floor(N/dse);
Xf = fft(in);

%zeroth order
tmp = ifft(Xf.*repmat(filts{1}.phi,1,L));
S0 = tmp(1:filts{1}.ds(J1+1):end,:);

%first layer
for j=1:J1
tmp = (ifft(Xf.*repmat(filts{1}.psi{j},1,L)));
scratch.X{1}(j,:,:) = tmp(1:filts{1}.ds(j):end,:);
end
U=abs(scratch.X{1});

UU=reshape_first(U,filts);
%now we have a K1 by K2 by Neff by L tensor
[K1,K2,Neff,L]=size(UU); 

%second layer: separable along 1st, 2nd and 3rd coordinates
%1st slice
Utmp=UU(:,:);
Uf = fft(Utmp);
Ubis = zeros(size(Utmp,1)*size(filts{2}.psi,2),size(Utmp,2));
rast=1;
for j=1:size(filts{2}.psi,2)
tmp = ifft( Uf.*repmat(filts{2}.psi{j},1,size(Uf,2)) );
tmp = tmp(1:filts{2}.ds(j):end,:);
st = size(tmp,1);
Ubis(rast:rast+st-1,:)=tmp;
scratch.rastbot(1,j)=rast;
scratch.rasttop(1,j)=rast+st-1;
rast=rast+st;
end
Ubis=Ubis(1:rast-1,:);
K1bis = size(Ubis,1);

%2nd slice
Ubis = reshape(Ubis,K1bis,K2,Neff,L);
Utmp = permute(Ubis,[2 1 3 4]);
Utmp = Utmp(:,:);
Uf = fft(Utmp);
Ubis = zeros(size(Utmp,1)*size(filts{3}.psi,2),size(Utmp,2));
rast=1;
for j=1:size(filts{3}.psi,2)
tmp = ifft( Uf.*repmat(filts{3}.psi{j},1,size(Uf,2)) );
tmp = tmp(1:filts{3}.ds(j):end,:);
st = size(tmp,1);
Ubis(rast:rast+st-1,:)=tmp;
scratch.rastbot(2,j)=rast;
scratch.rasttop(2,j)=rast+st-1;
rast=rast+st;
end
Ubis=Ubis(1:rast,:);
K2bis = size(Ubis,1);

%3rd slice
Ubis = reshape(Ubis,K2bis,K1bis,Neff,L);
Utmp = permute(Ubis,[3 2 1 4]);
Utmp = Utmp(:,:);
Uf = fft(Utmp);
Ubis = zeros(size(Utmp,1),size(filts{4}.psi,2),size(Utmp,2));
rast=1;
for j=1:size(filts{4}.psi,2)
tmp = ifft( Uf.*repmat(filts{4}.psi{j},1,size(Uf,2)) );
%tmp = tmp(1:filts{4}.ds(j):end,:);
Ubis(:,j,:)=tmp;
end
%scratch.X{2} = Ubis;

%phi output
tmp = ifft( Uf.*repmat(filts{4}.phi,1,size(Uf,2)));
%scratch.X{3} = tmp;
Utmp = abs(tmp);
Utmp=Utmp(1:filts{4}.dsphi:end,:);
Ulp = reshape(Utmp, size(Utmp,1),K1bis,K2bis,L);

%psi output
Utmp = abs(Ubis);
Utmp = Utmp(:,:);
Uf = fft(Utmp);
Utmp = ifft(Uf.*repmat(filts{4}.phi,1,size(Uf,2)));
Utmp = Utmp(1:filts{4}.dsphi:end,:);
Utmp = reshape(Utmp,size(Utmp,1),size(filts{4}.psi,2),K1bis,K2bis,L);
Utmp(:,size(filts{4}.psi,2)+1,:,:,:)=Ulp;
K3bis = size(filts{4}.psi,2)+1;

Utmp = reshape(Utmp, size(Utmp,1), K1bis*K2bis*K3bis, L);
%Utmp(:,end+1,:)=S0;
S = permute(Utmp,[2 1 3]);


