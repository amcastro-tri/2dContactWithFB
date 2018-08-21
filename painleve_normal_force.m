function lambda_N = painleve_normal_force(params, state)    
    m = params.m;
    half_length = params.half_length;
    g = params.g;
    mu = params.mu;
    
    % The column index is for solutions in time.
    psi = state(3, :);
    w = state(6, :);

    a = 1/m*(1+3*cos(psi).*(cos(psi)-mu*sin(psi)));
    b = g*(half_length*w.*w./g.*sin(psi)-1);

    lambda_N = zeros(size(state, 2), 1);    
    idx = a>0 & b < 0;
    lambda_N(idx) = -b(idx)./a(idx);
    
    %if(a>0 && b < 0)
    %    lambda_N = -b./a;
    %end
end