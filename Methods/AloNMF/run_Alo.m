javaaddpath tools
addpath (genpath('..'));
clear
rand('seed',1);

%%% set the size/sparsity of the data
m = 10000;
n = 8000;
k = 100;
rate = 0.3;

V = rand(n, m) * 100;
opts = getOptions();
W = rand(n, k) * max(max(V));
opts.H0 = rand(k, m);
opts.maxIter = 100;
opts.maxNumberThreads = 4;

[H, HIS] = AloProject(V, W, opts);



