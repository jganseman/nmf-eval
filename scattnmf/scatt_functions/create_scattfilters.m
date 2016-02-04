function filts = create_scattfilters(options)
%this function constructs the filters used for the time-frequency scattering
%transform. 

N=getoptions(options,'N',2^15);
T=getoptions(options,'T',2048);

%first layer filters
Q=getoptions(options,'Q',32);
os = getoptions(options,'os',1);
filt_opt.Q = Q;
filt_opt.J = T_to_J(T, filt_opt);
filt_opt.filter_format='fourier';
filt_opt.min_margin=0;
filt_opt.boundary='per';
filters=morlet_filter_bank_1d(N, filt_opt);
J=length(filters.psi.filter);
dse = getoptions(options,'dse',round(T*2^(-floor(J/Q)-os)));
Jinit=J;
for j= 1:J
filts{1}.psi{j}=filters.psi.filter{j};
filts{1}.ds(j) = dse;
end
filts{1}.phi = filters.phi.filter;
filts{1}.ds(J+1)=dse;% * 2^(floor(J/Q));
filts{1} = generate_dualfilters(filts{1});
filts{1}.Q=Q;
octave = filters.psi.meta.center/filters.psi.meta.center(1);
filts{1}.octaves = 1 - ceil(log2(octave));
rast=2;
filts{1}.tone(1)=1;
for r=2:J
if filts{1}.octaves(r) > filts{1}.octaves(r-1)
rast=1;
end
filts{1}.tone(r) = rast; rast=rast+1;
end
K1 = max(filts{1}.octaves);
K2 = filts{1}.Q;
for r=1:J
filts{1}.tonemap(r) = 1 + K1*(filts{1}.tone(r)-1)+filts{1}.octaves(r)-1;
end


%Haar filters along first level coefficients
Jhaar=getoptions(options,'Jhaar',3);
for j=1:Jhaar
slice=zeros(2^j,1);
slice(1:2^(j-1))=2^(-j);
slice(2^(j-1)+1:end)=-2^(-j);
filts{2}.haar{j}=slice;
filts{2}.ds(j)=max(1,2^(j-2));
end
filts{2}.haar{Jhaar+1}=2^(-Jhaar)*ones(2^Jhaar,1);
filts{2}.ds(Jhaar+1)=max(1,2^(Jhaar-2));

%second layer temporal filters
Neff = floor(N/dse);
J_2 = getoptions(options,'J_2',3);
clear filt_opt;
filt_opt.J=J_2;
filt_opt.filter_format='fourier';
filt_opt.min_margin=0;
filt_opt.boundary='per';
filters = morlet_filter_bank_1d(Neff, filt_opt);
J=length(filters.psi.filter);
for j=1:J
filts{4}.psi{j} = filters.psi.filter{j};
filts{4}.ds(j)=1;%2^(max(0,j-os));
end
filts{4}.phi = filters.phi.filter;
filts{4}.dsphi = 2^(J_2-1);% 2^(floor(Jinit/Q))
filts{4} = generate_dualfilters(filts{4});







