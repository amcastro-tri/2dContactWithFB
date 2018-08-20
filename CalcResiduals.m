function [rv, rl, Rvv, Rvl, Rlv, Rll] = CalcResiduals(v, lambda, problem_data, params)
M = problem_data.M;
Jn = problem_data.Jn;
Jt = problem_data.Jt;
pstar = problem_data.pstar;
h = params.h;
FB_lambda = params.FB_lambda;

vn = Jn * v;
vt = Jt * v;
[ft, dft_dvt, dft_dfn] = calc_friction_force(vt, lambda, params);

% Residual and Jacobians in v.
rv = M*v - pstar - h*Jn'*lambda - h*Jt'*ft;
Rvv =  M - h*Jt'*diag(dft_dvt)*Jt; % This is SPD! ==> invertible.
Rvl = -h * Jn' - h*Jt'*diag(dft_dfn);

% Residual and Jacobians in lambda.
eps = 1e-10;
rl = FischerBurmeisterLambda(vn, lambda, FB_lambda, eps);
Rlv = diag(FischerBurmeisterLambdaGradX(vn, lambda, FB_lambda, eps)) * Jn;
Rll = diag(FischerBurmeisterLambdaGradY(vn, lambda, FB_lambda, eps)); %Diagonal matrix. It might have zero diagonal elements.
 