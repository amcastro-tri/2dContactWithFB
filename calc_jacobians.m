function [Jn, Jt] = calc_jacobians(p_BoC_W)

nc = size(p_BoC_W, 2);

% Normal Jacobian.
Jn = zeros(nc, 3);
for ic = 1:nc
    Jn(ic, :) = [0.0, 1.0, p_BoC_W(1, ic)];
end

% Tangential Jacobian.
Jt = zeros(nc, 3);
for ic = 1:nc
    Jt(ic, :) = [1.0, 0.0, -p_BoC_W(2, ic)];
end


