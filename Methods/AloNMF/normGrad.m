function res = normGrad(gradX, X)
    res = norm(gradX(gradX < 0 | X > 0), 'fro')^2;
end