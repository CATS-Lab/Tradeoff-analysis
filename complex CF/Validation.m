syms kappa T t lambda

mu4_t= -(T-t)*exp(-sqrt(kappa)*(T-t))

mu1_t = simplify (1+kappa*int(mu4_t,t,T))

H = simplify(-(mu1_t+lambda*mu4_t))
LHS = diff(mu4_t,t) - (2*sqrt(kappa))*mu4_t-kappa*int(mu4_t,t,t,T);

RHS = 1;
result = simplify (LHS-RHS)