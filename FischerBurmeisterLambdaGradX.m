function z = FischerBurmeisterLambdaGradX(x, y, lambda)
    sx = x > 0;
    z = lambda * FischerBurmeisterGradX(x, y) - (1-lambda) * sx .* max(0, y);
end
