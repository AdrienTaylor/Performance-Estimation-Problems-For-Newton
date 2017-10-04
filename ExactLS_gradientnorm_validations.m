clear all; 
clc;

% syms kappa eps
kappa=sym('kappa','positive');
L=sym('L','positive');
eps=sym('eps','positive');
mu=kappa*L;
syms x0 x1 g0 g1;

%(1) interpolation
ineq1=1/(1-kappa)*((g0-g1)^2/L+mu*(x0-x1)^2-2*mu/L*(g0-g1)*(x0-x1))-...
    (g0-g1)*(x0-x1);

%(2) multipliers
lam=1+2*eps*kappa^(1/2)/((1-eps^2)^(1/2)*(1-kappa));
s11=3*eps/2-eps*(kappa+1/kappa)/4+(1-kappa)/2/...
    sqrt(kappa*(1-eps^2))-eps^2*(1-kappa)/sqrt(kappa*(1-eps^2));
s22=(eps*(1-kappa)-2*sqrt(kappa*(1-eps^2)))/...
    ((eps^2-1)*(1-kappa)-2*eps*sqrt(kappa*(1-eps^2)));
s21=eps*(1-kappa)/sqrt((kappa*(1-eps^2)))/2-1;

S=[s11 s21; s21 s22];

tau=(eps+(1-eps^2)^(1/2)*(1-kappa)/2/kappa^(1/2))^2;
%(3) exact line-search in inexact gradient direction
ineq2=g1*(x1-x0);
ineq3=-eps*g0^2*s11-eps*g1^2*s22-2*g1*g0*s21;

%(4) weighted sums: two formulations
expr=ineq1*lam*(L-mu)+ineq2*L*(1+kappa)+ineq3;
coeftot=kappa*(2*eps*kappa^(1/2)+(1-kappa)*(1-eps^2)^(1/2))...
    /((1-kappa)*(1-eps^2)^(1/2));
coefg1=eps*(1+kappa)...
    /(kappa^(1/2)*((1-eps^2)^(1/2)*(1-kappa)+2*eps*kappa^(1/2)));
coefg0=(1+kappa)/2/kappa;

expr_comp=g1^2-tau*g0^2+...
coeftot*(coefg1*g1-coefg0*g0+L*(x0-x1))^2;

%(5) verify the equivalence between the two formulations    
simplify((expr-expr_comp))% =0 so it is OK


