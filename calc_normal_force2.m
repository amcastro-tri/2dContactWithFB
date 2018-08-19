function [fn_p, dfdx, dfdxdot] = calc_normal_force2(x, xdot, k, d)

% Value of the normal force.
fn = k.*x + d.*xdot;
fn_p = max(0, fn);

% Gradient with respect to x, xdot.
Hfn = (fn >= 0);

dfdx = k.*Hfn;
dfdxdot = d.*Hfn;