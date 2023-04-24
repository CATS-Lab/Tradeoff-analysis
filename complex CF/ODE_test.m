clear
syms T t kappa lambda tau %theta phi 
theta = (tau*kappa+lambda)/2;
phi = sqrt(4*kappa-(tau*kappa+lambda)^2)/2;
%A = -i*(theta^2+phi^2)/(2*kappa*phi)*exp((theta-i*phi)*T);
%C = i*(theta^2+phi^2)/(2*kappa*phi)*exp((theta+i*phi)*T);
B = theta+i*phi;
D = theta-i*phi;
A = simplify (exp(-B*T)/(kappa*(1/D-1/B)));
C = simplify (-A*exp(B*T-D*T));
E=simplify(A/B);
F=simplify(C/D);
mu4_T = simplify(A*exp(B*T)+C*exp(D*T));
mu4_t = simplify(A*exp(B*t)+C*exp(D*t));

mu4_t1= simplify(-exp(-theta*(T-t))/phi*sin(phi*(T-t)));

%G = simplify(1+kappa*A/B*exp(B*T)+kappa*C/D*exp(D*T));
mu1_t = simplify(1+kappa*int(mu4_t, t, t,T));
alpha = atan(theta/phi);
%mu1_t1 = sqrt(kappa)*exp(-theta*(T-t))/phi*cos(phi*(T-t)-alpha);
mu1_t1 =exp(-theta*(T-t))/(2*phi)...
    *((phi+i*theta)*exp(-i*phi*(T-t))+(phi-i*theta)*exp(i*phi*(T-t)));
mu3_t = -mu1_t;

H = simplify (mu1_t+lambda*mu4_t);

beta = atan((theta-lambda)/phi)+(1-sign(theta-lambda))/2*pi;


%H1 = sqrt(kappa-tau*kappa*lambda)/phi*exp(-theta*(T-t))*cos(phi*(T-t)-beta);
H1 = sqrt(kappa-tau*kappa*lambda)/phi*exp(-theta*(T-t))*cos(phi*(T-t)-beta);

T_value =3;
t_value =1;
kappa_value =3;
lambda_value = 1;
tau_value =1;


mu4_t_value =vpa(subs(mu4_t,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]));
mu4_t1_value =vpa(subs(mu4_t1,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]));

mu4_t_diff = vpa(mu4_t_value-mu4_t1_value,2)

mu1_t_value =vpa(subs(mu1_t,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]));
mu1_t1_value =vpa(subs(mu1_t1,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]));

mu1_t_diff = vpa(mu1_t_value-mu1_t1_value,2)

H_value =vpa(subs(H,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]));

H1_value =vpa(subs(H1,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]));

H_diff = vpa(H_value-H1_value,2)

LHS = simplify(mu3_t+diff(mu4_t,t)-(tau*kappa+lambda)*mu4_t);
LHS_value =vpa(subs(LHS,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]))

LHS1 = simplify(mu1_t-1-kappa*int(mu4_t,t,t,T))
%LHS1_value =vpa(subs(LHS1,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]))

LHS2 = simplify(diff(mu1_t,t)+kappa*mu4_t)
LHS3 = simplify(diff(mu3_t,t)-kappa*mu4_t)

LH4 = -mu1_t-lambda*mu4_t;
RH4 = -1/(2*phi)*exp(-theta*(T-t))*((phi+i*(theta-lambda))*exp(-i*phi*(T-t))+(phi-i*(theta-lambda))*exp(i*phi*(T-t)));
Y4 = simplify(LH4-RH4)
Y4_value =vpa(subs(Y4,[T, t, kappa, lambda, tau],[T_value, t_value, kappa_value, lambda_value, tau_value]))
