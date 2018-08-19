function p_BoC_W = calc_contact_points(q, geometry)
%p_WBcm = q(1:2);
theta = q(3);

% The four corners, in the body frame.
p_BoC = geometry();
nc_max = size(p_BoC, 2);

% Rotation matrix
c = cos(theta);
s = sin(theta);
xhat = [c; s];
yhat = [-s; c];
R_WB = [xhat, yhat];

% The four corners, in the world frame.
p_BoC_W = zeros(2, nc_max);
for ic=1:nc_max
    p_BoC_W(:, ic) = R_WB * p_BoC(:, ic);
end
