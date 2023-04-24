clear
syms kappa tau T t t1 t2
theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));
H = sqrt(kappa)/phi*exp(-theta*(T-t))*sin(phi*(T-t)-alpha-pi/2);
mu2 = simplify(int(H,t,t1,t2))



