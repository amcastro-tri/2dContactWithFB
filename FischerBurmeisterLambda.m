function z = FischerBurmeisterLambda(x, y, lambda, eps)
  z = lambda * FischerBurmeister(x, y, eps) - (1-lambda) * max(0, x) .* max(0, y);
end
