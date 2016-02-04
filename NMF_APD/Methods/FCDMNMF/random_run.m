clear
max_iter = 500;
rand('seed',0);

%%% set the size/sparsity of the data
m = 500;
n = 1000;
k = 30;
rate = 0.3;

%%% generating V from random W and H
W_org = rand(n,k); W_org(rand(n,k)<rate)=0;
H_org = rand(k,m); H_org(rand(k,m)<rate)=0;
V = W_org * H_org + rand(n, m) * 10;

W = abs(rand(n,k));
H = abs(rand(k,m));
%% running GCD with trace=0, objGCD and timeGCD will just contain the final objective function and total running time
[W0 H0 objGCD timeGCD] = NMF_GCD(V,k,max_iter,W,H,0);

%% running GCD with trace=1, objGCD and timeGCD will contain a list of objective function and running time for each iteration
[W0 H0 objGCD timeGCD] = NMF_GCD(V,k,max_iter,W,H,1);


