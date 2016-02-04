function [gradW, gradH] = getGrad(A, W, H)
    gradW = W*(H*H') - A*H';
    gradH = (W'*W)*H - W'*A;
end
