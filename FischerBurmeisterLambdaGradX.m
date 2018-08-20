function z = FischerBurmeisterLambdaGradX(x, y, lambda, eps)
    sx = x > 0;
    z = lambda * FischerBurmeisterGradX(x, y, eps) - (1-lambda) * sx .* max(0, y);
end
