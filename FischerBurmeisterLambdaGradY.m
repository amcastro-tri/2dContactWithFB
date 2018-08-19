function z = FischerBurmeisterLambdaGradY(x, y, lambda)
    sy = y > 0;
    z = lambda * FischerBurmeisterGradY(x, y) - (1-lambda) * sy .* max(0, x);
end
