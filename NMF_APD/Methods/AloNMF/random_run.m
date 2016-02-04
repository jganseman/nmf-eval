addpath (genpath('..'));
clear
rand('seed',1);

%%% set the size/sparsity of the data
m = 10000;
n = 500;
k = 30;
rate = 0.3;

V = rand(n, m) * 100;

opts = getOptions();
opts.W0 = rand(n, k) * max(max(V));
opts.H0 = rand(k, m);
opts.maxIter = 200;
maxNumCompThreads(4)
%% running GCD with trace=0, objGCD and timeGCD will just contain the final objective function and total running time
[W, H, HIS] = AlNMF(V, opts);
norm(V -W*H, 'fro')^2/2.0;

%% running GCD with trace=1, objGCD and timeGCD will contain a list of objective function and running time for each iteration
%[W, H, HIS] = AlNMF(V, k, opts);


