syms s kappa lambda tau 
eq1 = exp(s.*tau)+(lambda+kappa.*tau)*s+kappa == 0;
sol = solve(eq1,s);
sol2 =subs(sol,[kappa,lambda],[0.01,0.01]);
taus =[0.001:1:50];
sol3 = subs(sol2,tau,taus);
sol4 = vpa(simplify(sol3));
figure(1);clf;hold on
plot(taus,sol4);
figure(2);clf;hold on

plot(xs,lambertw(0,xs));


