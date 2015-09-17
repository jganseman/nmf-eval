function [dout] = audioscatt_bkw_haar(dSin, dUin, scratch, filts, options)
%here we compute the backprop of the two-layer scattering operator for speech separation

dSin = dSin .* scratch.PSaux ;

dSin = permute(dSin, [2 1 3]);

dU = reshape(dSin,scratch.Utmpsize1);
%we split the lowpass from the rest of the channels
dUlp = squeeze(dU(:,end,:,:));
dU = dU(:,1:end-1,:,:);

dU = reshape(dU, scratch.Utmpsize2);
dUtmp=zeros(scratch.Utmpsize3);
dUtmp(1:filts{4}.dsphi:end,:) = dU;

dU = ifft(fft(dUtmp).*repmat(conj(filts{4}.phi),1,size(dUtmp,2)));
dU = reshape(dU, scratch.Utmpsize4);
dU = dU .* scratch.PShp;

dUlp = reshape(dUlp, scratch.Utmpsize5);
dUtmp = zeros(scratch.Utmpsize3(1),scratch.Utmpsize5(2));
dUtmp(1:filts{4}.dsphi:end,:) = dUlp ;
dUtmp = dUtmp .* scratch.PSlp;
dUtmp = ifft(fft(dUtmp).*repmat(conj(filts{4}.phi),1,size(Uf,2)));

dUf = zeros(size(dU,1),size(dU,3));

for j=1:size(filts{4}.psi,2)
tmp = squeeze(dU(:,j,:));
dUf = dUf + ifft(fft(tmp).*repmat(conj(filts{4}.psi{j},1,size(dUf,2))));
end

%join phi and psi outputs
dUf = dUf + dUtmp ;
dUf = reshape(dUf, scratch.Utmpsize6);
dUf = permute(dUf, [2 1 3]);

if scratch.dohaar
aux = zeros(scratch.Utmpsize7);
dUf = reshape(dUf, scratch.Ubissize);
for j=1:scratch.Jhaar
slice = dUf(scratch.rastbot(1,j):scratch.rasttop(1,j),:);
tmp=zeros(scratch.M0h{j},size(slice,2));
tmp(1:filts{2}.ds(j):end,:) = slice;
aux = aux + conv2(tmp,filts{2}.chaar{j},'full');
end
dU = reshape(aux, scratch.Usize);
else
dU = reshape(dUf, scratch.Usize);
end


%%%%
dU = dU + dUin;

dU = dU .* scratch.P ;







