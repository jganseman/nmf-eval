function [W0, H0] = initWH(V, K)
    maxV = max(max(V));
    [N, M] = size(V);
    W0 = rand(N, K) .* maxV;
    H0 = rand(K, M);
end