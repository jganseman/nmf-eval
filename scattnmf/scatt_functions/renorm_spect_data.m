function [out,norms]=renorm_spect_data(in, stds, epsilon)
out=in;

out = out./repmat(stds,1,size(out,2));
norms = sqrt(sum(abs(out).^2)) + epsilon;
out = out./repmat(norms,size(out,1),1);


