function [state_x, fn_all, ft_all, vn_all, vt_all, vn_err, vt_err] = box_discrete_update(itime, state_x0, params)

% Params
m = params.m;
I = params.I;
g = params.g;
h = params.h;
stiction_tolerance = params.stiction_tolerance;
relative_tolerance = params.relative_tolerance;
ev = relative_tolerance * stiction_tolerance;

% State
q0 = state_x0(1:3);
v0 = state_x0(4:6);
p_WBo = q0(1:2);
%theta = q(3);
%v_WBo = v(1:2);
%w = v(3);

% Sizes
nv = 3;

p_BoC_W = calc_contact_points(q0, params.geometry);
nc_max = size(p_BoC_W, 2);

% Calc signed distance at t = t0.
x0_all = zeros(nc_max,1);
for ic = 1:nc_max
    p_WC = p_WBo + p_BoC_W(:, ic);
    x0_all(ic) = -p_WC(2);
end
idx = find(x0_all > -1e-4);
x0 = x0_all(idx);
nc = length(x0);
p_BoC_W = p_BoC_W(:, idx);

[Jn, Jt] = calc_jacobians(p_BoC_W);

% Genralized forces
tau = [0; -m*g; 0]; % + Jn'*fn + Jt'*ft;


% Newton-Rapshon loop
max_iters = 100;
v = v0;
lambda = zeros(nc, 1);
M = [m, 0, 0;
     0, m, 0;
     0, 0, I];
pstar = M*v0 + h*tau;

% Problem data
problem_data.M = M;
problem_data.Jn = Jn;
problem_data.Jt = Jt;
problem_data.pstar = pstar;

% Relaxploation seems to help.
w = 0.9;

vn_err = 2*ev;
vt_err = 2*ev;
for it=1:max_iters 
    % Normal/tangential velocities
    vn = Jn*v;
    vt = Jt*v;
    
    % Residuals and Jacobians.
    [rv, rl, Rvv, Rvl, Rlv, Rll] = CalcResiduals(v, lambda, problem_data, params);
    
    % Factorize and then solve in Matlab?
    % We can use Cholesky for Rvv.
    Rvvi = inv(Rvv);
    
    Rvvi_rv = Rvvi*rv;
    A = Rll - Rlv*Rvvi*Rvl;  
    p = -rl + Rlv*Rvvi_rv;
    
    % lambda might be overconstrained. We use LSQ solution.
    dl = A\p; %How to call explicitly in Matlab?
    
    % v should be unique (physics).
    dv = -Rvvi_rv - Rvvi*Rvl*dl;           
          
    dvn = Jn * dv;
    dvt = Jt * dv;
    
    % Limit velocity update.
    alpha = 1.0;
    for ic=1:nc
        vt0 = vt(ic);
        
        % We started in the "strong" gradients region. Safe.
        if (abs(vt0) < stiction_tolerance)
            continue
        end
        
        vt1 = vt(ic) + dvt(ic);
        
        % If velocity changes direction.
        if (vt0 * vt1 < 0)
            % Clamp to within the stiction region, keeping the sign.
            vt_alpha = sign(vt0) * stiction_tolerance/2;
            
            alpha_ic = (vt_alpha-vt0)/dvt(ic);
            if (alpha_ic < 0 || alpha_ic > 1)
                error('Alpha is outsude [0, 1]');
            end
            
            alpha = min(alpha, alpha_ic);
        end        
    end
    
    % Update velocities and labmdas
    v = v + w*alpha*dv;
    lambda = lambda + w*alpha*dl;
    
    % Normal/tangential velocities
    vn = Jn*v;
    vt = Jt*v;
    
    % Errors
    vn_err = norm(dvn);
    vt_err = norm(dvt);
    
    % Check if converged in velocities ONLY (since lambda might still
    % change for overconstrained systems even if bounded by the generalized
    % contact forces).
    if (vn_err < ev && vt_err < ev) 
        % Compute fores for reporting by LSQ.
        tau_n = Jn' * lambda;
        fn = pinv(Jn')*tau_n;
        %fn = lambda;
        ft = calc_friction_force(vt, fn, params);
        break;
    end        
end

if (vn_err > ev || vt_err > ev)
    % If we are here is because the NR iteration failed. Abort.
    msg = sprintf('NR iteration did not converge.\n It: %d.\n vn_err: %g.\n vt_err: %g. \n Time step: %d\n', it, vn_err, vt_err, itime);
    error(msg);
end

% Update the state
q = q0 + h*v;
state_x = [q; v];

fn_all = zeros(nc_max, 1);
ft_all = zeros(nc_max, 1);
vn_all = zeros(nc_max, 1);
vt_all = zeros(nc_max, 1);

fn_all(idx) = fn;
ft_all(idx) = ft;
vn_all(idx) = vn;
vt_all(idx) = vt;


