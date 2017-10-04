%% Validation code, convergence in distance to optimality
clc; clear all;

syms x0 xs g0 d f0 fs
syms gam mu L eps R

gs=0;
x1=x0-gam*d;

% all inequalities are in the form (...)<=0

%(1) interpolation
ineq1=f0-fs+g0*(xs-x0)+...
    1/(2*(1-mu/L))*(1/L*(g0-gs)^2+mu*(x0-xs)^2-2*mu/L*(x0-xs)*(g0-gs));
ineq2=fs-f0+gs*(x0-xs)+...
    1/(2*(1-mu/L))*(1/L*(g0-gs)^2+mu*(x0-xs)^2-2*mu/L*(x0-xs)*(g0-gs));

%(2) inexact search direction
ineq3=(d-g0)^2-eps^2*(g0)^2;

% % initial condition (not needed for proof)
% ineq4=(g0)^2-R^2;

%(3) multipliers 
rho=(1-(1-eps)*mu*gam);
tau=rho^2;

lam0=2*rho*gam*(1-eps);
lam1=gam*rho/mu/eps;

%(4) weighted sums: two formulations
expr_tot=lam0*(ineq1+ineq2)+lam1*(ineq3);

coefsos1=gam*mu^2*(1-eps)*(2-gam*(1-eps)*(L+mu))/(L-mu);
coefsos1d=(L-mu)/((1-eps)*mu^2*(2-gam*(1-eps)*(L+mu)));
coefsos1g0=-(L+mu)*(1-gam*mu*(1-eps))/(mu^2*(2-gam*(1-eps)*(L+mu)));
sos1=coefsos1*(x0-xs+coefsos1d*d+coefsos1g0*g0)^2;

coefsos2=gam*(1-gam*mu*(1-eps))*(2*mu-eps*(L+mu)-gam*mu*(1-eps)*(L+mu))...
    /(eps*mu^2*(1-eps)*(2-gam*(1-eps)*(L+mu)));
coefsos2g0=-(1-eps);
sos2=coefsos2*(d+coefsos2g0*g0)^2;
expr_reworked=(x1-xs)^2-tau*(x0-xs)^2+sos1+sos2;

%(5) verify the equivalence between the two formulations    
simplify(expr_tot-expr_reworked)%=0: OK !
