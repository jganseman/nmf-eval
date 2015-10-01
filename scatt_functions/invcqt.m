function out=invcqt(in, options, filts)

os = getoptions(options,'os',1);
Q = getoptions(options,'Q',32);

[N,M, L]=size(in);
%1) upscale modulus 
J=size(filts.psi,2);
dse = filts.ds(1);%getoptions(options,'dse',round(T*2^(-floor(J/Q)-os)));
M0 = M * dse;

out=zeros(M0,L);

[xx0,yy0]=meshgrid(1:dse:M0,1:N);
[xx,yy]=meshgrid(1:M0,1:N);

for l=1:L
in_up(:,:,l)=interp2(xx0,yy0,in(:,:,l),xx,yy,'spline');
end

iters = getoptions(options, 'maxiters', 20);

%init with random phase (except low frequencies)
z = randn(size(in_up)) + i * randn(size(in_up));
z(:,end,:) = 0;

for it=1:iters

fprintf('it %d, dist is %f \n', it, norm(abs(z(:)) - in_up(:)))
z = in_up.*exp(i*angle(z));
out=0*out;
%3) pseudoinverse
for l=1:L

coeffs = z(:,:,l);
coeffs=transpose(coeffs);
for j=1:J
out(:,l) = out(:,l) + real(ifft(fft(coeffs(:,j)).*filts.dpsi{j}));
end
out(:,l) = out(:,l) + real(ifft(fft(coeffs(:,J+1)).*filts.dphi));

fout = fft(out);
for j=1:J
coeffs(:,j) = ifft(fout.*filts.psi{j});
end
coeffs(:,J+1) = 0*ifft(fout.*filts.phi);
z(:,:,l) = transpose(coeffs);
end
end



