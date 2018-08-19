function z = FischerBurmeisterLambda(x, y, lambda)
  z = lambda * FischerBurmeister(x, y) - (1-lambda) * max(0, x) .* max(0, y);
end
