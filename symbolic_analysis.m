clear;

syms delta kappa tau a_bar

theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));

A = sym(zeros (2,4));
B = sym(zeros (1,4));
A_dot = sym(zeros (2,4));
B_dot = sym(zeros (1,4));

denom = exp(theta*delta).*cos(phi*delta)-1;
numer = exp(theta*delta).*sin(phi*delta);
beta = asin (numer/sqrt(denom^2+numer^2));

deltas = [delta,pi/phi-delta,delta,pi/phi-delta];
b = [0.5* a_bar,0,-0.5* a_bar,0];
c = [0,delta* a_bar,delta* a_bar,0];
 

for j = 1: 4
    t = deltas(j);
    A(1,j) = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
    A(2,j) = exp(-theta*t)*sin(phi*t)*1/phi;  
    Bb = 2*tau*t-1/theta*tau*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
     -1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
    Bc = -tau*(exp(-theta*t)*cos(phi*t)-1)- theta*tau/phi*exp(-theta*t)*sin(phi*t);
    B(j)= Bb*b(j)+Bc*c(j);
    
    A_dot(1,j) = -kappa/phi*exp(-theta*t)*sin(phi*t);
    A_dot(2,j)  = -1/phi*theta*exp(-theta*t)*sin(phi*t)+exp(-theta*t)*cos(phi*t);
    Bb_dot = -2*tau*(exp(-theta*t)*cos(phi*t)-1)+1/phi*(2-2*theta*tau)*exp(-theta*t)*sin(phi*t);
    Bc_dot = tau*kappa/phi*exp(-theta*t)*sin(phi*t);
    B_dot(j)= Bb_dot*b(j)+Bc_dot*c(j);
end
D = zeros(1,2);
D_dot = zeros(1,2);


E = B(3);
E_dot = B_dot(3);
for l3 = [1,2]
    D_temp3 = A(l3,3);
    D_dot_temp3 = A_dot(l3,3);    
    E_temp3 = A(l3,3);
    E_dot_temp3 = A_dot(l3,3);
    
    if l3 == 1
        E = E+E_temp3*B(2);
        E_dot = E_dot+E_dot_temp3*B(2);
    else
        E = E+E_temp3*B_dot(2);
        E_dot = E_dot+E_dot_temp3*B_dot(2);
    end
    for l2 = [1,2]
        if l3 == 1           
            D_temp2 = D_temp3*A(l2,2);
            D_dot_temp2 = D_dot_temp3*A(l2,2);            
            E_temp2 = E_temp3*A(l2,2);
            E_dot_temp2 = E_dot_temp3*A(l2,2);
        else
            D_temp2 = D_temp3 *A_dot(l2,2);
            D_dot_temp2 = D_dot_temp3 *A_dot(l2,2);
            E_temp2 = E_temp3*A_dot(l2,2);
            E_dot_temp2 = E_dot_temp3*A_dot(l2,2);
        end
        if l2 == 1
            E = E+E_temp2*B(1);
            E_dot = E_dot+E_dot_temp2*B(1);
        else
            E = E+E_temp2*B_dot(1);
            E_dot = E_dot+E_dot_temp2*B_dot(1);
        end
        for l1 = [1,2]
            if l2 == 1
                D_temp1 = D_temp2 *A(l1,1);
                D_dot_temp1 = D_dot_temp2 *A(l1,1);
                E_temp1 = E_temp2 *A(l1,1);
                E_dot_temp1 = E_dot_temp2 *A(l1,1);
            else
                D_temp1 = D_temp2 *A_dot(l1,1);
                D_dot_temp1 = D_dot_temp2 *A_dot(l1,1);
                E_temp1 = E_temp2 *A_dot(l1,1);
                E_dot_temp1 = E_dot_temp2 *A_dot(l1,1);
            end
            
            if l1 ==1
                D = D+ D_temp1 *[A(1,4),A(2,4)];
                D_dot = D_dot+ D_dot_temp1 *[A(1,4),A(2,4)];
                E = E+E_temp1*B(4);
                E_dot = E_dot+E_dot_temp1*B(4);
            else
                D = D+ D_temp1 *[A_dot(1,4),A_dot(2,4)];
                D_dot = D_dot+ D_dot_temp1 *[A_dot(1,4),A_dot(2,4)];
                E = E+E_temp1*B_dot(4);
                E_dot = E_dot+E_dot_temp1*B_dot(4);
            end
        end
    end
end

% ks = [1:K];
% t3 = -pi/phi + delta + ks*2*pi/phi;
% Delta3 = E;
% Delta_dot_3 = E_dot;
% for k = 2:K
%     Delta3(k) = D(1)*Delta3(k-1)+D(2)*Delta_dot_3(k-1)+E;
%     Delta_dot_3(k) = D_dot(1)*Delta3(k-1)+D_dot(2)*Delta_dot_3(k-1)+E_dot;
% end



Coef_1 =  1/(1-D(1))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
Coef_2 = D(2)/((1-D(1))*(1-D_dot(2)))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
Delta3_inf = Coef_1*E + Coef_2*E_dot;

Coef_dot_1 = D_dot(1)/((1-D(1))*(1-D_dot(2)))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
Coef_dot_2 = 1/(1-D_dot(2))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));

Delta3_dot_inf = Coef_dot_1*E + Coef_dot_2*E_dot;

t = (2*alpha+beta)/phi - delta;
A1 = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
A2 = exp(-theta*t)*sin(phi*t)*1/phi;  
Bb = 2*tau*t-1/theta*tau*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
 -1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
Bc = -tau*(exp(-theta*t)*cos(phi*t)-1)- theta*tau/phi*exp(-theta*t)*sin(phi*t);
B4= Bb*b(4)+Bc*c(4);

Delta4_inf =A1*Delta3_inf + A2*Delta3_dot_inf+B4;


D = simplify(D)
D_dot = simplify(D_dot)
E = simplify(E)
%Delta3_inf = simplify(Delta3_inf)
%Delta4_inf = simplify(Delta4_inf)