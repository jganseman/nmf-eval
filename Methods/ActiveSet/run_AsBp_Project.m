clear
rand('seed',1);
addpath ..
%%% set the size/sparsity of the data
m = 10000;
n = 8000;
k = 100;
rate = 0.3;

V = rand(n, m) * 1000;
%opts = getOptions();
W = rand(n, k) * max(max(V));
opts.H0 = rand(k, m);
opts.maxIter = 100;
opts.maxNumberThreads = 4;

[H, HIS] = ActiveSetProject(V, W, opts);

[H, HIS] = BlockpivotProject(V, W, opts);



