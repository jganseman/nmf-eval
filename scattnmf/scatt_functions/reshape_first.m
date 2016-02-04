function UU=reshape_first(U, filts)
%this function maps the cqt coefficients into a 2D array to perform the second layer transform

[J, N, L] = size(U);

%filts{1}.octave and filts{1}.tone;
K1 = max(filts{1}.octaves);
K2 = filts{1}.Q;

U= U(:,:);

Utmp = zeros(K1*K2, N*L);

Utmp(filts{1}.tonemap,:)=U;

UU=reshape(Utmp,K1,K2,N,L);






