function [W, iter] = doIter(HH, hh, W0, maxIter, tolerence, verbose)
    [k, ~] = size(HH);
    [~, n] = size(hh);
    sqrtH = diag(HH);
    sqrtH(sqrtH==0)=1;
    for i=1:k,
        HH(:,i) = HH(:,i)./(sqrtH*sqrtH(i));
    end
    iter = 0;
    for i=1:n,
        hh(:,i) = hh(:,i)./sqrtH;
        W0(:,i) = W0(:,i)./sqrtH;
        [W0(:,i), it] = optimize(HH, hh(:,i), W0(:,i), maxIter, tolerence, verbose);
        iter = iter + it;
        W0(:,i) = W0(:,i).*sqrtH;
    end
    W = W0;
end

function [x, iter] = optimize(Q, q, x0, maxIter, tolerence, verbose)
    x = x0;
    df = Q*x + q;
    [k, ~] = size(Q);
    dbarf = df;
    dbarf(df > 0 & x == 0) = 0;
    preDerivative = norm(dbarf, 'fro')^2;
    for iter=1:2,
        for i=1:k,
            dbarf = df;
            dbarf(df > 0 & x == 0) = 0;
            [~, pos] = max(dbarf);
            dx = max(0, x(pos) - df(pos)) - x(pos);
            df = df + dx * Q(:,pos);
            x(pos) = x(pos) + dx;
            %fprintf('%E\n', x'*Q*x/2 + q'*x);
        end
%         dbarf = df;
%         dbarf(df >= 0 & x == 0) = 0;
%         den = dbarf'*Q*dbarf;
%         if (den <= 1e-100), break; end
%         alpha = dbarf'*dbarf/den;
%         if (alpha ~= alpha || alpha > 1e50), break; end
%         x = max(0, x - alpha * dbarf);
%         df = Q*x + q;
        dbarf(df > 0 & x == 0) = 0;
        derivative = norm(dbarf, 'fro')^2;
        %fprintf('%E\n', x'*Q*x/2 + q'*x);
        if (derivative <= preDerivative * tolerence), break; end
    end
    %fprintf('finished\n');
end