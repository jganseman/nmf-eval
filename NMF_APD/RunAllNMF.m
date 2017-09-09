addpath (genpath('Methods'));
javaaddpath java;

datasets = {'faces', 'digits', 'cifar-10'};
Ks = [60 80 100];

datasets = {'digits'};
Ks = [8];

rand('seed', 1);

% for i=1:length(datasets),
%     X = MatrixReader.readMatrix(sprintf('datasets/%s/train.gz', datasets{i}))';
%     [N, M] = size(X);
%     K = Ks(i);
%     opts = getOptions();
%     [opts.W0, opts.H0] = initWH(X, K);
%     writeMatrix(sprintf('datasets/%s/W0_%d.csv', datasets{i}, K), opts.W0', ',');
%     writeMatrix(sprintf('datasets/%s/H0_%d.csv', datasets{i}, K), opts.H0, ',');
% end
%@QnNMF, @NtNMF, 
maxNumCompThreads(4)
Methods = {@LeeNMF, @PGLIN, @HALSacc, @NMF_GCD, @ActiveSet, @Blockpivot, @NeNMF, @AlNMF}
%@PGLIN
Methods = {@AloExactNMF, @Blockpivot}
Methods = {@AlNMF}
%Methods = {@Blockpivot}

opts = getOptions();
opts.maxIter = 30;
opts.tolerance = 1e-2;
opts.maxNumberThreads = maxNumCompThreads;
    
for i=1:length(datasets),
    fprintf('%s\n', datasets{i});
    X = MatrixReader.readMatrix(sprintf('datasets/%s/train.gz', datasets{i}))';
    [N, M] = size(X);
    K = Ks(i);
    opts.W0 = csvread(sprintf('datasets/%s/W0_%d.csv', datasets{i}, K))';
    opts.H0 =csvread(sprintf('datasets/%s/H0_%d.csv', datasets{i}, K));
    for m=1:length(Methods)
        fprintf('%s\n', func2str(Methods{m}));
        [W, H, HIS] = Methods{m}(X, opts);
        s = sprintf('results/%s_%s_%d', datasets{i}, func2str(Methods{m}), K);
        writeMatrix(sprintf('%s_W.cvs', s), W, ',');
        writeMatrix(sprintf('%s_H.cvs', s), H, ',');
        %size((0:opts.maxIter)'),
        writeMatrix(sprintf('%s_log.cvs', s), [(0:opts.maxIter)' HIS.time' HIS.iterW' HIS.iterH' HIS.error' HIS.gradW' HIS.gradH' HIS.grad'], ',');
        fprintf('%s    %.10f\n', func2str(Methods{m}), norm(X-W*H, 'fro')^2/2);
        %[H, HIS] = AloProject(X, W, opts);
    end
end 
