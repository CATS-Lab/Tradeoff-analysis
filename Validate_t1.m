clear
syms T t kappa tau delta

theta = (tau*kappa)/2;

phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));

H =sqrt(kappa)/phi*exp(-theta*(T-t))*cos(phi*(T-t)-(alpha+pi));

varphi = i*phi;
mu4_t = 1/(2*varphi)*(exp(-(theta+varphi)*(T-t))-exp(-(theta-varphi)*(T-t)));
mu1_t = simplify(1+kappa*int(mu4_t, t, t,T));
H1 = simplify (-mu1_t);

H_diff = H1 - H;

denom = exp(theta*delta)*cos(phi*delta)-1;
numer = exp(theta*delta)*sin(phi*delta);
beta = acot(denom/numer);


t1 = T -(2*alpha+beta)/phi;

check = int(H,t,t1,t1+delta);

digits(64)
tau_value = vpa(0.5);
kappa_value = 0.9*4/tau_value^2;
T_value = 50;
t_value = 10;
phi_value = vpa(subs(phi,[kappa, tau],[kappa_value, tau_value]));
theta_value = vpa(subs(theta,[kappa, tau],[kappa_value, tau_value]));
delta_value = 0.3295; %1.7* pi/phi_value;


denom_value = vpa(subs(denom,[T, kappa, tau, delta],[T_value, kappa_value, tau_value,delta_value]));
numer_value = vpa(subs(numer,[T, kappa, tau, delta],[T_value, kappa_value, tau_value,delta_value]))

beta_value = vpa(subs(beta,[T, kappa, tau, delta],[T_value, kappa_value, tau_value,delta_value]))





H_diff_value = vpa(subs(H_diff,[T, kappa, tau, t],[T_value, kappa_value, tau_value,t_value]));
check_value = vpa(subs(check,[T, kappa, tau, delta],[T_value, kappa_value, tau_value,delta_value]))


t1_value = vpa(subs(t1,[T, kappa, tau, delta],[T_value, kappa_value, tau_value,delta_value]));
t1_value = vpa(mod(t1_value,pi/phi_value))
