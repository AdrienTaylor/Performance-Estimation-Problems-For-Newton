%% Validation code, convergence in objective function accuracy
clc; clear all;

syms f0 f1 fs g0 g1 x0 xs d
syms rho gam eps L mu
gs=0;
x1=x0-gam*d;

% all inequalities are in the form (...)<=0

%(1) interpolation
ineq1=f1-f0+g1*(x0-x1)+1/2/L*(g0-g1)^2+mu/(2*(1-mu/L))...
    *(x0-x1-1/L*(g0-g1))^2;
ineq2=f0-fs+g0*(xs-x0)+1/2/L*(g0-gs)^2+mu/(2*(1-mu/L))...
    *(x0-xs-1/L*(g0-gs))^2;
ineq3=f1-fs+g1*(xs-x1)+1/2/L*(g1-gs)^2+mu/(2*(1-mu/L))...
    *(x1-xs-1/L*(g1-gs))^2;

%(2) inexact search direction
ineq4=(d-g0)^2-eps^2*(g0)^2;

%(3) multipliers 
lam01=rho;
lams0=rho*(1-rho);%\lambda_{*0}
lams1=1-rho;%\lambda_{*1}
lam2=rho*gam/2/eps;
tau=rho^2;

%(4) weighted sums: two formulations
expr=ineq1*lam01+ineq2*lams0+ineq3*lams1+ineq4*lam2;


coefsos1=(gam*rho*(L*(-2*eps*gam*mu+rho-1)+mu*(rho+1)))...
    /(2*eps*(L*(rho-1)+mu*(rho+1)));
coefsos1g0=((eps+1)*L*(rho-1)-(eps-1)*mu*(rho+1))...
    /(L*(2*eps*gam*mu-rho+1)-mu*(rho+1));

coefsos2=L*mu*(1-rho^2)/2/(L-mu);
coefsos2d=-gam/(rho+1);
coefsos2g0=-rho/(mu*rho+mu);
coefsos2g1=-1/(mu*rho+mu);

coefsos3=(rho*(L+mu)-(L-mu))/(2*mu*(rho+1)*(L-mu));
coefsos3d=2*gam*L*mu*rho/(L*(rho-1)+mu*(rho+1));
coefsos3g0=rho*((L*(rho-1)-mu*(rho+1)))/(L*(rho-1)+mu*(rho+1));

coefsos4=rho*((eps+1)*gam*L-1-rho)*((eps-1)*gam*mu+1-rho)...
    /(L*(2*eps*gam*mu+1-rho)-mu*(rho+1));

expr_tot=f1-fs-tau*(f0-fs)+coefsos1*(d+g0*coefsos1g0)^2+...
    coefsos2*(x0-xs+d*coefsos2d+g0*coefsos2g0+g1*coefsos2g1)^2+...
    coefsos3*(g1+d*coefsos3d+g0*coefsos3g0)^2-...
    coefsos4*(g0)^2;
    
%(5) verify the equivalence between the two formulations    
simplify(expr-expr_tot)



