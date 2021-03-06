%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Painleve
rod_length = 0.2;
m = 0.1;
g = 9.81;
I = m * rod_length * rod_length / 12.0;
stiction_tolerance = 1.0e-4;
relative_tolerance = 1e-6; %INVESTIGATE!: Some convergence issues observed for small time step (1e-5) and tole=1e-3 at the singularity in fn.
FB_lambda = 1.0;
mu = 1.5;  % > 4/3
theta0 = -20 / 180 * pi;
w0 = -1.5 * (2 * pi);
vx0 = 2.0;
vy0 = -rod_length/2*w0*cos(-theta0); % Notice w0 has the sign of phi (negative)
y0 = rod_length/2*sin(-theta0);
sim_time = 0.6;
h = 5.0e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m = 0.1;
% g = 9.81;
lengths = [0.05, 0.1];
% I = m * (lengths(1)^2 + lengths(2)^2)/12.0;
% stiction_tolerance = 1.0e-5;
% relative_tolerance = 1e-3;
% FB_lambda = 0.9;
% mu = 0.35;
% y0 = 0.3;
% vx0 = -1.5;
% vy0 = 0.0;
% w0 = 0.0;
% sim_time = 2.0;
% h = 1.0e-3;

% Define the geometry for a box.
nc_max = 4;
p_BoC = zeros(2, nc_max);
p_BoC(:, 1) = [-lengths(1); -lengths(2)] / 2;
p_BoC(:, 2) = [ lengths(1); -lengths(2)] / 2;
p_BoC(:, 3) = [ lengths(1);  lengths(2)] / 2;
p_BoC(:, 4) = [-lengths(1);  lengths(2)] / 2;
box = @() p_BoC;

%Define the geometry for a box with multicontact (overconstrained).
nc_max = 8;
p_BoC = zeros(2, nc_max);
p_BoC(:, 1) = [-lengths(1); -lengths(2)] / 2;
p_BoC(:, 2) = [0; -lengths(2)] / 2;
p_BoC(:, 3) = [ lengths(1); -lengths(2)] / 2;
p_BoC(:, 4) = [ lengths(1); 0] / 2;
p_BoC(:, 5) = [ lengths(1);  lengths(2)] / 2;
p_BoC(:, 6) = [0;  lengths(2)] / 2;
p_BoC(:, 7) = [-lengths(1);  lengths(2)] / 2;
p_BoC(:, 8) = [-lengths(1);  0] / 2;
box2 = @() p_BoC;

% Define the geometry for an n-edges polygon.
num_sides = 2;
radius = 0.05;
t = (1/num_sides/2:1/num_sides:1)*2*pi;
nc_max = num_sides;
p_BoC = zeros(2, nc_max);
p_BoC(1, :) = cos(t) * radius;
p_BoC(2, :) = sin(t) * radius;
polygon = @() p_BoC;

% Define the geometry for a rod
nc_max = 2;
p_BoC = zeros(2, nc_max);
p_BoC(:, 1) = [-rod_length; 0]/2;
p_BoC(:, 2) = [rod_length; 0]/2;
rod = @() p_BoC;


% Save parameters into a struct
params.m = m;
params.I = I;
params.g = g;
params.stiction_tolerance = stiction_tolerance;
params.relative_tolerance = relative_tolerance;
params.FB_lambda = FB_lambda;
params.mu = mu;
params.h = h;
params.geometry = rod;

x0 = [0; y0; theta0; 
      vx0; vy0; w0];

nsteps = ceil(sim_time/h);
xx = zeros(nsteps, 6);
fn = zeros(nsteps, nc_max);
ft = zeros(nsteps, nc_max);
tt = zeros(nsteps, 1);
vn = zeros(nsteps, nc_max);
vt = zeros(nsteps, nc_max);
xp = zeros(nsteps, nc_max);
vn_err = zeros(nsteps, 1);
vt_err = zeros(nsteps, 1);
for it=1:nsteps
    tt(it) = it * h;
    [x, fn(it,:), ft(it,:), vn(it,:), vt(it,:), vn_err(it), vt_err(it)] = box_discrete_update(it, x0, params);
    
%     if (any(fn(it,:) > 0))
%         fn(it,:)
%         ft(it,:)
%         
%         error('Contact at step: %d\n', it);
%     end
    
    xx(it, :) = x;
    x0 = x;
end
