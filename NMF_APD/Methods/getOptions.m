function opts = getOptions()
    opts.maxIter = 300;
    opts.verbose = 1;
    opts.type = 'PLAIN';
    opts.tolerance = 1e-4;
    opts.timeLimit = 1e10;
    opts.alpha = 0.01;
    opts.delta = 0.01;
    opts.maxNumberThreads = 4;
end