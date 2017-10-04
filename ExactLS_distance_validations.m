clear all; 
clc;

% syms kappa eps
kappa=sym('kappa','positive');
L=sym('L','positive');
eps=sym('eps','positive');
mu=kappa*L;
syms x0 x1 xs g0 g1;
gs=0;

% all inequalities are in the form (...)<=0

%(1) interpolation

ineq1=1/(1-kappa)...
    *((g0-gs)^2/L+mu*(x0-xs)^2-2*kappa*(g0-gs)*(x0-xs))+g0*(xs-x0);
ineq2=1/(1-kappa)...
    *((g1-gs)^2/L+mu*(x1-xs)^2-2*kappa*(g1-gs)*(x1-xs))+g1*(xs-x1);

%(2) multipliers

lam0=(1-kappa)/mu*(1-2*eps^2+...
    eps*(1-eps^2)^(1/2)/2/kappa^(1/2)/(1-kappa)*(-1-kappa^2+6*kappa));
lam1=1/mu-1/L;
lam2=1/mu+1/L;
s11=(3*eps/2-eps*(kappa+1/kappa)/4+(1-kappa)/2/sqrt(kappa*(1-eps^2))-...
    eps^2*(1-kappa)/sqrt(kappa*(1-eps^2)))/(L*mu);
s22=((eps*(1-kappa)-2*sqrt(kappa*(1-eps^2)))...
    /((eps^2-1)*(1-kappa)-2*eps*sqrt(kappa*(1-eps^2))))/(L*mu);
s21=(eps*(1-kappa)/sqrt((kappa*(1-eps^2)))/2-1)/(L*mu);

tau=(eps+(1-eps^2)^(1/2)*(1-kappa)/2/kappa^(1/2))^2;
S=[s11 s21; s21 s22];

%(3) exact line-search in inexact gradient direction
ineq3=g1*(x1-x0);
ineq4=-eps*g0^2*s11-eps*g1^2*s22-2*g1*g0*s21;

%(4) weighted sums: two formulations
expr=ineq1*lam0+ineq2*lam1+ineq3*lam2+ineq4;

coeftot=(2*eps*(1-eps^2)^(1/2)*kappa^(1/2)+(1-eps^2)*(1-kappa))...
    /(kappa*(1-kappa));
coefx0=-(1+kappa)/2;
coefg0=(1-eps*(1-kappa)/(2*(1-eps^2)^(1/2)*kappa^(1/2)))/L;
coefg1=(1-kappa)...
    /(2*eps*(1-eps^2)^(1/2)*kappa^(1/2)+(1-eps^2)*(1-kappa))/L;



expr_comp=(x1-xs)^2-tau*(x0-xs)^2+...
    coeftot*(coefx0*(x0-xs)+coefg0*g0+coefg1*g1)^2;

%(5) verify the equivalence between the two formulations    
simplify((expr-expr_comp))% =0 so it is OK


