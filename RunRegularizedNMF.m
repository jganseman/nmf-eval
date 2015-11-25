addpath (genpath('Methods'));
javaaddpath java;

datasets = {'faces', 'digits', 'cifar-10'};
Ks = [60 80 100];

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
Methods = {@LeeNMF, @PGLIN, @Blockpivot, @ActiveSet, @NMF_GCD, @NeNMF,  @HALSacc, @AlNMF}
%@PGLIN
%@NeNMF, @ActiveSet,
Methods = {@NeNMF, @ActiveSet, @AlNMF};
Methods = {@AlNMF};
%Methods = {@Blockpivot}

for i=1:length(datasets),
    fprintf('%s\n', datasets{i});
    X = MatrixReader.readMatrix(sprintf('datasets/%s/train.gz', datasets{i}))';
    [N, M] = size(X);
    K = Ks(i);
    opts = getOptions();
    opts.maxIter = 30;
    %opts.mu_1 = 0.001;
    opts.mu_2 = 0.001;
    %opts.beta_1 = 0.001;
    opts.beta_2 = 0.001;
    opts.maxNumberThreads = maxNumCompThreads;
    opts.W0 = csvread(sprintf('datasets/%s/W0_%d.csv', datasets{i}, K))';
    opts.H0 =csvread(sprintf('datasets/%s/H0_%d.csv', datasets{i}, K));
    for m=1:length(Methods)
        fprintf('%s\n', func2str(Methods{m}));
        [W, H, HIS] = Methods{m}(X, opts);
        s = sprintf('results/%s_%s_%d', datasets{i}, func2str(Methods{m}), K);
        writeMatrix(sprintf('%s_W_regularized_2.cvs', s), W, ',');
        writeMatrix(sprintf('%s_H_regularized_2.cvs', s), H, ',');
        %size((0:opts.maxIter)'),
        writeMatrix(sprintf('%s_log_regularized_2.cvs', s), [(0:opts.maxIter)' HIS.time' HIS.iterW' HIS.iterH' HIS.error' HIS.gradW' HIS.gradH' HIS.grad'], ',');
        fprintf('%s    %.10f\n', func2str(Methods{m}), norm(X-W*H, 'fro')^2/2);
        %[H, HIS] = AloProject(X, W, opts);
        if strcmp(func2str(Methods{m}), func2str(@AlNMF)) == 0,
            continue;
        end;
        opts = rmfield (opts, 'mu_2');
        %writeMatrix(sprintf('%s_W_regularized_mu_2.cvs', s), W, ',');
        %writeMatrix(sprintf('%s_H_regularized_mu_2.cvs', s), H, ',');
        %size((0:opts.maxIter)'),
        writeMatrix(sprintf('%s_log_regularized_mu_2.cvs', s), [(0:opts.maxIter)' HIS.time' HIS.iterW' HIS.iterH' HIS.error' HIS.gradW' HIS.gradH' HIS.grad'], ',');
        fprintf('%s    %.10f\n', func2str(Methods{m}), norm(X-W*H, 'fro')^2/2);
        %[H, HIS] = AloProject(X, W, opts);
    end
end 
