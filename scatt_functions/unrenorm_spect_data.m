function out=unrenorm_spect_data(in, stds,norms)

out = in.*repmat(norms,size(in,1),1);
out = out.*repmat(stds,1,size(in,2));


