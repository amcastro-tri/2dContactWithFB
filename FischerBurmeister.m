function z = FischerBurmeister(x, y, eps)
  e2 = eps^2;
  z = sqrt(x.^2 + y.^2+e2) - x - y;
end


