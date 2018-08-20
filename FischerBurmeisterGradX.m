function z = FischerBurmeisterGradX(x, y, eps)
    %eps2 = 1e-10;
    %if (norm([x; y]) < eps2)
    %    z = -1;
    %else
    %    z = x ./ (sqrt(x.^2 + y.^2)) - 1;
    %end
    e2 = eps^2;
    z = x ./ (sqrt(x.^2 + y.^2+e2)) - 1;    
end
