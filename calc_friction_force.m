function [ft, dft_dvt, dft_dfn] = calc_friction_force(vt, fn, params)

nc = length(fn);

mu = params.mu;
stiction_tolerance = params.stiction_tolerance;
relative_tolerance = params.relative_tolerance;
ev = relative_tolerance * stiction_tolerance;
ev2 = ev*ev;

ft = zeros(nc, 1);
dft_dvt = zeros(nc, 1);
dft_dfn = zeros(nc, 1);
for ic=1:nc
    vt_ic = vt(ic);
    slip = abs(vt_ic);
    mu_ic = stribeck_friction2(slip, mu, stiction_tolerance);
    
    sign = vt_ic / sqrt(vt_ic^2 + ev2);
    
    ft(ic) = -mu_ic * fn(ic) * sign;
    
    % Calc dft/dfn
    dft_dfn(ic) = -mu_ic * sign;
    
    % Calc dft/dvt
    mu_prime = stribeck_friction_prime2(slip, mu, stiction_tolerance); %dmu/dvt
    dft_dvt(ic) = -mu_prime*fn(ic);    
end