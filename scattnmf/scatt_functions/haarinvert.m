function out=haarinvert(U, filts,scratch)

[K, L] = size(U);

bord=4;
J=size(filts.haar,2);

for j=1:J
slice = U(scratch.rastbot(1,j):scratch.rasttop(1,j),:);
%upscale to original resolution
dse = filts.ds(j);
Mtmp=size(slice,1);
M0= scratch.M0h{j}; %Mtmp * dse; 
[xx0,yy0]=meshgrid(1:L,1:dse:M0);
%add border
xx0bis=zeros(2*bord+size(xx0,1),size(xx0,2));
yy0bis=zeros(2*bord+size(xx0,1),size(xx0,2));
xx0bis(1:bord,:)=repmat(xx0(1,:),bord,1);
xx0bis(end-bord+1:end,:)=repmat(xx0(end,:),bord,1);
xx0bis(bord+1:end-bord,:)=xx0;
rr=dse*[1:bord]';
yy0bis(1:bord,:)=repmat(yy0(1,1)-flipud(rr),1,size(yy0bis,2));
yy0bis(end-bord+1:end,:)=repmat(yy0(end,1)+rr,1,size(yy0bis,2));
yy0bis(bord+1:end-bord,:)=yy0;
slicebis=0*xx0bis;
slicebis(1:bord,:)=repmat(slice(1,:),bord,1);
slicebis(end-bord+1:end,:)=repmat(slice(end,:),bord,1);
slicebis(bord+1:end-bord,:)=slice;
[xx,yy]=meshgrid(1:L,1:M0);
chunk{j}=interp2(xx0bis,yy0bis,slicebis,xx,yy,'cubic',0);
end

%this is OK
out=haarinv(chunk);




