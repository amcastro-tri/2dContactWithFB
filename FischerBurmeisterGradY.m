function z = FischerBurmeisterGradY(x, y)
    z = y ./ (sqrt(x.^2 + y.^2)+1e-14) - 1;
end
