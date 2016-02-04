%test filters are OK


x = randn(Npad,1);

J1 = size(filts{1}.psi,2);
xf = fft(x);

reco = 0*x;
for j=1:J1
   transx(:,j) =  ifft(xf.*filts{1}.psi{j});
   reco = reco + ifft(fft(transx(:,j)).*filts{1}.dpsi{j});
end
transx(:,J1+1) = ifft(xf.*filts{1}.phi);
reco = real( reco + ifft(fft(transx(:,J1+1)).*filts{1}.dphi));

