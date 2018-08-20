function z = FischerBurmeisterLambdaGradY(x, y, lambda, eps)
    sy = y > 0;
    z = lambda * FischerBurmeisterGradY(x, y, eps) - (1-lambda) * sy .* max(0, x);
end
