
close all; clear; clc;

syms t psi(t) theta(t) phi(t) 

L1h = [cos(psi(t)) sin(psi(t)) 0 ; -sin(psi(t)) cos(psi(t)) 0 ; 0 0 1];

L21 = [cos(theta(t)) 0 sin(theta(t)) ; 0 1 0 ; -sin(theta(t)) 0 cos(theta(t))];

Lb2 = [1 0 0 ; 0 cos(phi(t)) -sin(phi(t)) ; 0 sin(phi(t)) cos(phi(t))];

Lbh = Lb2*L21*L1h;

Lhb = transpose(Lbh);

Lhb_p = diff(Lhb,t,1);

RES = simplify(Lbh*Lhb_p);

display(RES);

%%

syms t u(t) I t0 t0p s1 s2
delta_m = s1*t + s2*t^2;
theta = 1/I *(s1*t^3/6 + s2*t^4/12) + t0p*t + t0;
w = ((delta_m - theta) - diff(u,t,1))/diff(theta,t,1);
u_sol = dsolve(delta_m^2 - theta^2 == 2*(diff(w,t,1) - u*diff(theta,t,1)));
pretty(u_sol)








