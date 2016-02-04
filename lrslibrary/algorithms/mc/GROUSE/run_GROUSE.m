
[numr,numc] = size(M);
I = randi([0 1],numr,numc); % ones(size(M));
maxrank = 1;
maxCycles = 100;
step_size = 0.1;

[Usg, Vsg, err_reg] = grouse(M,I,numr,numc,maxrank,step_size,maxCycles);
L = Usg*Vsg';
S = M - L;

% show_2dvideo(M,m,n);
% show_2dvideo(M.*I,m,n);
% show_2dvideo(L,m,n);
% show_2dvideo(S,m,n);
