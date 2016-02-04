%% [struct] = process_matrix(string, string, 2dmatrix, struct = []);
%
function results = process_matrix(method_id, algorithm_id, M, opts)
  method_name = get_method_name_by_id(method_id);
  algorithm_name = get_algorithm_name_by_id(algorithm_id);
  displog(['Running ' method_name ' with ' algorithm_name]);
  
  %%% RPCA methods
  % i.e: results = process_matrix('RPCA', 'FPCP', M, []);
  if(strcmp(method_id,'RPCA'))
    results = run_algorithm_rpca(algorithm_id, M, opts);
  end
  
  %%% ST methods
  % i.e: results = process_matrix('ST', 'GRASTA', M, []);
  if(strcmp(method_id,'ST'))
    results = run_algorithm_st(algorithm_id, M, opts);
  end
  
  %%% MC methods
  % i.e: results = process_matrix('MC', 'GROUSE', M, []);
  if(strcmp(method_id,'MC'))
    results = run_algorithm_mc(algorithm_id, M, opts);
  end
  
  %%% LRR methods
  % i.e: results = process_matrix('LRR', 'FastLADMAP', M, []);
  if(strcmp(method_id,'LRR'))
    results = run_algorithm_lrr(algorithm_id, M, opts);
  end
  
  %%% TTD methods
  % i.e: results = process_matrix('TTD', '3WD', M, []);
  if(strcmp(method_id,'TTD'))
    results = run_algorithm_ttd(algorithm_id, M, opts);
  end
  
  %%% NMF methods
  % i.e: results = process_matrix('NMF', 'NMF-MM', M, []);
  if(strcmp(method_id,'NMF'))
    results = run_algorithm_nmf(algorithm_id, M, opts);
  end
  
  %%% Apply hard thresholding
  results.O = hard_threshold(results.S); % imagesc(results.O);
  
  displog('Decomposition finished');
end
