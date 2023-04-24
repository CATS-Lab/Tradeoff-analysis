clear
syms T t kappa lambda tau va
theta = (tau*kappa+lambda)/2;
varphi = sqrt((tau*kappa+lambda)^2-4*kappa)/2;
%A = -i*(theta^2+varphi^2)/(2*kappa*varphi)*exp((theta-i*varphi)*T);
%C = i*(theta^2+varphi^2)/(2*kappa*varphi)*exp((theta+i*varphi)*T);

mu4_t = 1/(2*varphi)*(exp(-(theta+varphi)*(T-t))-exp(-(theta-varphi)*(T-t)));
mu1_t = simplify(1+kappa*int(mu4_t, t, t,T));
H = simplify (-mu1_t-lambda*mu4_t);

A = (theta-lambda-varphi)*(theta-varphi)*(exp((theta+varphi)*va)-1);
B = (theta-lambda+varphi)*(theta+varphi)*(exp((theta-varphi)*va)-1);
t11 = T-1/(2*varphi)*log(A/B);


LHS = int(H,t,t11,t11+va);
%LHS = simplify(diff(mu4_t,t) - (kappa*tau+lambda)*mu4_t-kappa*int(mu4_t,t,t,T));

T_value =12;
t_value =1;
tau_value = 0.9;
lambda_value = 0.9*(1/tau_value);
kappa_value = 0.9*(2-lambda_value*tau_value-2*sqrt(1-lambda_value*tau_value))/tau_value^2;


va_value =2;
theta_value = vpa(subs(theta,[T, t, kappa, lambda, tau,va],[T_value, t_value, kappa_value, lambda_value, tau_value,va_value]))
varphi_value = vpa(subs(varphi,[T, t, kappa, lambda, tau,va],[T_value, t_value, kappa_value, lambda_value, tau_value,va_value]))

check_sqrt = (kappa_value*tau_value+lambda_value)^2-4*kappa_value

check = vpa(theta_value-lambda_value+varphi_value)

t11_value =vpa(subs(t11,[T, t, kappa, lambda, tau,va],[T_value, t_value, kappa_value, lambda_value, tau_value,va_value]))

LHS_value =vpa(subs(LHS,[T, t, kappa, lambda, tau,va],[T_value, t_value, kappa_value, lambda_value, tau_value,va_value]))

