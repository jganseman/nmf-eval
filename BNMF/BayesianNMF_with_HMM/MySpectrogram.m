function [Magnitude,Phase,Power]=MySpectrogram(x, alen, ulen)
%Spectrogram of an input waveform
% x is the speech vector
% alen is the analysis frame length, ulen is the update length
%
% Naser Mohammadiha November 2010


x=[zeros(alen-ulen,1);x];
naf=floor(size(x,1)/ulen);  % Number of analysis frames
x=[x;zeros(alen-ulen,1)];
x=x(1:naf*ulen+alen-ulen);

n1 = 1;n2 = alen;
M=alen/2+1;

Magnitude=zeros(M,naf);
Power=zeros(M,naf);
Phase=zeros(M,naf);
win=hann(alen,'periodic');
win=win/sqrt(sum(win.^2)); %normalize window

for n=1:naf % Counter over analysis frames
    xf = x(n1:n2);
    xf=xf.*win;
    X=fft(xf);
    sft=abs(X(1:M));
    
    
    Magnitude(:,n)=sft;
    Power(:,n)=Magnitude(:,n).^2;
    Phase(:,n)=angle(X(1:M));
    n1 = n1 + ulen;
    n2 = n2 + ulen;
end
end