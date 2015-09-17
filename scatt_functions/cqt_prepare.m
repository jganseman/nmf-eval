function [filts, filt_opt]= cqt_prepare(options)

N=getoptions(options,'N',4096);
T=getoptions(options,'T',2048);
Q=getoptions(options,'Q',32);

filt_opt.Q = Q;
filt_opt.J = T_to_J(T, filt_opt);
filt_opt.filter_format='fourier';
filt_opt.min_margin=0;
filt_opt.boundary='per';
filters=morlet_filter_bank_1d(N, filt_opt);

for j= 1:length(filters.psi.filter)
filts.psi{j}=filters.psi.filter{j};
end
filts.phi = filters.phi.filter;
filts = generate_dualfilters(filts);

filts.Q=Q;
filts.N=N;




