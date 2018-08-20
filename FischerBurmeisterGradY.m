function z = FischerBurmeisterGradY(x, y, eps)   
    %eps2 = 1e-10;
    %if (norm([x; y]) < eps2)
    %    z = -1;
    %else
    %    z = y ./ (sqrt(x.^2 + y.^2)) - 1;
    %end
    e2 = eps^2;
    z = y ./ (sqrt(x.^2 + y.^2 + e2)) - 1;
end
