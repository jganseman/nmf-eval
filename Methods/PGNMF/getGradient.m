function [gradW,gradH] = getGradient(A,W,H,type,alpha,beta)
    switch type
        case 'plain'
            gradW = W*(H*H') - A*H';
            gradH = (W'*W)*H - W'*A;
        case 'regularized'
            gradW = W*(H*H') - A*H' + alpha*W;
            gradH = (W'*W)*H - W'*A + beta*H;
        case 'sparse'
            k=size(W,2);
            betaI = beta*ones(k,k);
            gradW = W*(H*H') - A*H' + alpha*W;
            gradH = (W'*W)*H - W'*A + betaI*H;
    end
end