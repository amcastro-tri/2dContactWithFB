function z = FischerBurmeisterGradX(x, y)
    z = x ./ (sqrt(x.^2 + y.^2)+1e-14) - 1;
end
