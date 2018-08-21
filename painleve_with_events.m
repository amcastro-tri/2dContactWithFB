function [tt, xx, fn, ft, aa, bb, te,ye,ie, params] = painleve_with_events

% Rod's properties.
rod_length = 0.2;
half_length = rod_length/2;
m = 0.1;
g = 9.81;
I = m * rod_length * rod_length / 12.0;
mu = 1.5;  % > 4/3

h = 1e-4; %desired interval for output.

params.half_length = half_length;
params.m = m;
params.g = g;
params.mu = mu;
params.h = h;

% Initial conditions (satisfying contraint vn=0 at t=0).
theta0 = 20 / 180 * pi;
w0 = 1.5 * (2 * pi);
vx0 = 2.0;
vy0 = rod_length/2*w0*cos(theta0);
y0 = rod_length/2*sin(theta0);

% Time stepping.
sim_time = 0.25;

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@events);

x0 = [0; y0; theta0; vx0; vy0; w0]';
[t1,y1,te1,ye1,ie1] = ode45(@f,0:h:sim_time,x0,options);

% We stopped at an event. Now compute impact event.
y_plus = calc_impact(ye1);

% Keep integratin until sim_tim.
[t2,y2,te2,ye2,ie2] = ode45(@f,te1:h:sim_time,y_plus,options);
%t2=[]; y2=[]; te2=[]; ye2=[]; ie2=[];

% Output into format ready to animate.
tt = [t1; t2];
xx = [y1; y2];

[lambda, aa, bb] = normal_force(xx');

xx(:, 3) = -xx(:, 3);

te = [te1; te2];
ye = [ye1; ye2];
ie = [ie1; ie2];

fn = [zeros(length(tt), 1) lambda];
ft = -mu * fn;

% Define the geometry for the rod
nc_max = 2;
p_BoC = zeros(2, nc_max);
p_BoC(:, 1) = [-rod_length; 0]/2;
p_BoC(:, 2) = [rod_length; 0]/2;
rod = @() p_BoC;

params.geometry = rod;

% -----------------------------------------------------------------------
% Nested functions -- problem parameters provided by the outer function.
%

function state = calc_impact(state0)
    xx0 = state0(1);
    yy0 = state0(2);
    psi0 = state0(3);
    vvx0 = state0(4);
    vvy0 = state0(5);
    ww0 = state0(6);

    % Angular velocity can be computed right from conservation of angular
    % moement about the impact point. This is EQUIVALENT to STICTION AT
    % IMPACT.
    L0 = I*ww0 + m*half_length*sin(psi0)*vvx0 + m*half_length*cos(psi0)*vvy0;    
    I_ell = I + m*half_length*half_length;

    % Conservation of angular momentum.
    w = L0/I_ell;       
    
    vx = half_length*w*sin(psi0);  % stiction at impact.
    vy = half_length*w*cos(psi0);  % elastic collision (e=0).

    lambdaT = m*(vx-vvx0);
    lambdaN = m*(vy-vvy0);

    lambdaT
    mu*lambdaN

    E0 = 0.5*m*(vvx0^2+vvy0^2) + 0.5*I*ww0^2
    E = 0.5*m*(vx^2+vy^2) + 0.5*I*w^2
    
    state = [xx0 yy0 psi0 vx vy w]';
end

function [lambda_N, a, b] = normal_force(state)
    %nt = size(state, 1);
    %state = repmat(state,           
    
    % The column index is for solutions in time.
    y = state(2, :);
    psi = state(3, :);
    w = state(6, :);

    gy = y - half_length*sin(psi);

    a = 1/m*(1+3*cos(psi).*(cos(psi)-mu*sin(psi)));
    b = g*(half_length*w.*w./g.*sin(psi)-1);

    lambda_N = zeros(size(state, 2), 1);    
    % The last condition is used just to avoid forces during free flight.
    idx = a>0 & b < 0 & gy < 1e-4;
    lambda_N(idx) = -b(idx)./a(idx);
    
    %if(a>0 && b < 0)
    %    lambda_N = -b./a;
    %end
end

function dydt = f(t,state)
    xc = state(1);
    yc = state(2);
    psi = state(3);
    vx = state(4);
    vy = state(5);
    w = state(6);

    lambda_N = normal_force(state);
    lambda_T = -mu * lambda_N;

    ax = lambda_T/m;
    ay = lambda_N/m - g;
    aw = (-lambda_T*half_length*sin(psi) - lambda_N*half_length*cos(psi))/I;

    dydt = [vx vy w ax ay aw]';
end

% -----------------------------------------------------------------------

function [value,isterminal,direction] = events(t,state)
    psi = state(3);
    w = state(6);

    a = 1/m*(1+3*cos(psi)*(cos(psi)-mu*sin(psi)));
    b = g*(half_length*w*w/g*sin(psi)-1);

    value = a;
    isterminal = 1;
    direction  = -1;
end

% -----------------------------------------------------------------------

end  % painleve_with_events