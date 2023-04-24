D1 = 0.35;
D2 = 0.2;
D1_dot = 0.4;
D2_dot = 0.45;

E = 13;
E_dot = 15;

Delta = 0;
Delta_dot = 0;

Deltas = 0;
Deltas_dot = 0;

max_iter = 50;
for i =1: max_iter
    Delta = D1*Delta+D2*Delta_dot+E;
    Delta_dot = D1_dot *Delta + D2_dot*Delta_dot +E_dot;
    Deltas = [Deltas, Delta];
    Deltas_dot = [Deltas_dot, Delta_dot];
end

Coef_1 =  1/(1-D1)*1/(1-D1_dot*D2/((1-D1)*(1-D2_dot)));
Coef_2 = D2/((1-D1)*(1-D2_dot))*1/(1-D1_dot*D2/((1-D1)*(1-D2_dot)));
Delta_inf = Coef_1*E + Coef_2*E_dot;

Coef_dot_1 = D1_dot/((1-D1)*(1-D2_dot))*1/(1-D1_dot*D2/((1-D1)*(1-D2_dot)));
Coef_dot_2 = 1/(1-D2_dot)*1/(1-D1_dot*D2/((1-D1)*(1-D2_dot)));

Delta_dot_inf = Coef_dot_1*E + Coef_dot_2*E_dot;

figure(1);cla;hold all;
plot([0:max_iter],Deltas, 'r','linewidth',2);
plot([0,max_iter],[Delta_inf,Delta_inf], '--r');
plot([0:max_iter],Deltas_dot, 'b','linewidth',2);
plot([0,max_iter],[Delta_dot_inf,Delta_dot_inf], '--b');



