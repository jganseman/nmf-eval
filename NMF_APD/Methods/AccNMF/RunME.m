% Run me : Example on the cbcl face dataset
clc; load cbclim; 

% Generating random initial iterates
[m,n] = size(M); r = 49; 
U0 = rand(m,r); V0 = rand(r,n);
maxiter = 300; timelimit = 5000000;

disp('************************************************************************');
disp('Comparaison of MU, HALS and PG algorithms with their accelerated variants');
disp('on the CBCL face dataset with factorization rank 49.');
disp('************************************************************************');

% Original MU
%[Um,Vm,em,tm] = MUacc(M,U0,V0,0,0,maxiter,timelimit); 
%disp(sprintf('original MU terminated with final error %f',em(end))); 
% Accelerated MU
%[Uma,Vma,ema,tma] = MUacc(M,U0,V0,2,0.1,maxiter,timelimit);
%disp(sprintf('accelerated MU terminated with final error %f',ema(end))); 
% Original HALS
%[Uh,Vh,eh,th] = HALSacc(M,U0,V0,0,0,maxiter,timelimit);
%disp(sprintf('original HALS terminated with final error %f',eh(end))); 
% Accelerated HALS
%[Uha,Vha,eha,tha] = HALSacc(M,U0,V0,0.5,0.1,maxiter,timelimit);
%disp(sprintf('accelerated HALS terminated with final error %f',eha(end))); 
% Projected gradient of Lin
[Up,Vp,ep,tp] =  PGLIN(M,U0,V0,1e-16,timelimit,maxiter);    
disp(sprintf('original PG terminated with final error %f',ep(end))); 
% Accelerated projected gradient of Lin
[Upa,Vpa,epa,tpa] = PGLINacc(M,U0,V0,0.5,0,maxiter,timelimit);
disp(sprintf('accelerated PG terminated with final error %f',epa(end)));        

% Plot evolution of the error w.r.t. time
figure; plot(tm,em,'b'); hold on; plot(tma,ema,'b--'); 
plot(th,eh,'k');  plot(tha,eha,'k--'); 
plot(tp,ep,'r');  plot(tpa,epa,'r--'); 
xlabel('time (s.)'); ylabel('||M-UV||_F');
title('blue: MU, black: HALS, red: PG, plain: original, dashed: accelerated'); 

% Display of basis elements obtained with different algorithms
lig = 7; Li = 19; Co = 19;
affichage(M(:,1:50:end),lig,Li,Co); title('Sample of images from CBCL dataset'); 
affichage(Um,lig,Li,Co); title('MU basis'); 
affichage(Uma,lig,Li,Co); title('accelerated MU basis'); 
affichage(Uh,lig,Li,Co); title('HALS basis'); 
affichage(Uha,lig,Li,Co); title('accelerated HALS basis'); 
affichage(Up,lig,Li,Co); title('PG basis'); 
affichage(Upa,lig,Li,Co); title('accelerated PG basis'); 