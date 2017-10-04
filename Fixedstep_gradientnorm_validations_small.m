%% Validation code, convergence in gradient norm
clc; clear all;

syms x0 g0 g1 d f0 f1
syms gam mu L eps R

x1=x0-gam*d;

% all inequalities are in the form (...)<=0

%(1) interpolation
ineq1=f0-f1+g0*(x1-x0)+...
    1/(2*(1-mu/L))*(1/L*(g0-g1)^2+mu*(x0-x1)^2-2*mu/L*(x0-x1)*(g0-g1));
ineq2=f1-f0+g1*(x0-x1)+...
    1/(2*(1-mu/L))*(1/L*(g0-g1)^2+mu*(x0-x1)^2-2*mu/L*(x0-x1)*(g0-g1));

%(2) inexact search direction
ineq3=(d-g0)^2-eps^2*(g0)^2;

% %initial condition (not really needed for proof)
% ineq4=(g0)^2-R^2;

%(3) multipliers 
rho=(1-(1-eps)*mu*gam);
tau=rho^2;

lam0=2*rho/gam/(1-eps);
lam1=gam*mu*rho/eps;

%(4) weighted sums: two formulations
expr_tot=lam0*(ineq1+ineq2)+lam1*(ineq3);

coefNorm1=((2-(1-eps)*gam*(L+mu))/((1-eps)*gam*(L-mu)));
coefNorm1d=(gam*(L+mu)*((eps-1)*gam*mu+1))/((eps-1)*gam*(L+mu)+2);
coefNorm1g0=-2*((eps-1)*gam*mu+1)/((eps-1)*gam*(L+mu)+2);
coefNorm1g1=1;

coefNorm2=(1-eps)*gam*(1-(1-eps)*gam*mu)*(-(1-eps)*gam*mu*(L+mu)-...
    eps*(L+mu)+2*mu)/(eps*(2-(1-eps)*gam*(L+mu)));
coefNorm2d=1/(eps-1);
coefNorm2g0=1;
coefNorm2g1=0;

expr_reworked=(g1)^2-tau*(g0)^2+...
    coefNorm1*(d*coefNorm1d+g0*coefNorm1g0+g1*coefNorm1g1)^2+...
    coefNorm2*(d*coefNorm2d+g0*coefNorm2g0+g1*coefNorm2g1)^2;

%(5) verify the equivalence between the two formulations    
simplify(expr_reworked-expr_tot) % OK (=0)
