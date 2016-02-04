function [Ufinal] = audioreconstruct2(S, Plp, Php, filts, scratch)

U= permute(S,[2 1 3]);
[Nfin, rien, L]=size(U);
%U= reshape(U,Nfin,K2bis,K1bis,L);

%start with the psi outputs
%lp is the last coordinate
%1) upscale 
Ubis = U(:,:);
dse = filts{4}.dsphi;
M0 = scratch.M0;
N0 = size(Ubis,2);
K1bis = scratch.K1bis;
K2bis = scratch.K2bis;

[xx0,yy0]=meshgrid(1:N0,1:dse:M0);
[xx,yy]=meshgrid(1:N0,1:M0);
Ubis=interp2(xx0,yy0,Ubis,xx,yy,'spline');

U=reshape(Ubis,M0,K2bis,K1bis,L);

Uhp = U(:,1:K2bis-1,:,:);
Ulp = squeeze(U(:,K2bis,:,:));

Uf = fft(Uhp(:,:));
Uhp = ifft(Uf.*repmat(filts{4}.dphi,1,size(Uf,2)));
Uhp = reshape(Uhp, M0, K2bis-1,K1bis,L);

%introduce the recorded phases
Uhp = Uhp .* Php;
Ulp = Ulp .* Plp;

%invert the temporal wavelets
Ubis = zeros(M0,K1bis,L);
for j=1:K2bis-1
tempo = Uhp(:,j, :, :);
tf = fft(tempo(:,:));
Ubis = Ubis + real(ifft( tf .* repmat(filts{4}.dpsi{j},1,size(tf,2))));
end
tf = fft(Ulp(:,:));
Ubis = Ubis + real(ifft( tf .* repmat(filts{4}.dphi,1,size(tf,2))));


%invert the Haar transform
Uh = permute(Ubis,[2 1 3]);

[K1bis, M0, L] = size(Uh);

if scratch.dohaar
Ubis = haarinvert(Uh(:,:), filts{2}, scratch);
else
Ubis = Uh;
end

Ufinal = reshape(Ubis,size(Ubis,1), M0, L);







