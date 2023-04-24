clear;

syms d k t a

theta = (t*k)/2;
phi = sqrt(4*k-(t*k)^2)/2;
alpha = asin(theta/sqrt(k));
d = pi/phi;


A = sym(zeros (2,4));
B = sym(zeros (1,4));
A_dot = sym(zeros (2,4));
B_dot = sym(zeros (1,4));

 denom = exp(theta*d).*cos(phi*d)-1;
 numer = exp(theta*d).*sin(phi*d);
beta = asin (numer/sqrt(denom^2+numer^2))+pi;

%beta = atan(phi/theta);

ds = [d,pi/phi-d,d,pi/phi-d];
b = [0.5* a,0,-0.5* a,0];
c = [0,d* a,d* a,0];
 

for j = 1: 4
    interval = ds(j);
    A(1,j) = exp(-theta*interval)*cos(phi*interval) + theta/phi*exp(-theta*interval)*sin(phi*interval);
    A(2,j) = exp(-theta*interval)*sin(phi*interval)*1/phi;  
    Bb = -(2*t*interval-1/theta*t*(1-k*t^2)*(exp(-theta*interval)*cos(phi*interval)-1)...
     -1/phi*t*(3-k*t^2)*exp(-theta*interval)*sin(phi*interval));
    Bc = -(-t*(exp(-theta*interval)*cos(phi*interval)-1)- theta*t/phi*exp(-theta*interval)*sin(phi*interval));
    B(j)= Bb*b(j)+Bc*c(j);
    
    A_dot(1,j) = -k/phi*exp(-theta*interval)*sin(phi*interval);
    A_dot(2,j)  = -1/phi*theta*exp(-theta*interval)*sin(phi*interval)+exp(-theta*interval)*cos(phi*interval);
    Bb_dot = -(-2*t*(exp(-theta*interval)*cos(phi*interval)-1)+1/phi*(2-2*theta*t)*exp(-theta*interval)*sin(phi*interval));
    Bc_dot = -(t*k/phi*exp(-theta*interval)*sin(phi*interval));
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




Coef_1 =  1/(1-D(1))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
Coef_2 = D(2)/((1-D(1))*(1-D_dot(2)))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
d3_inf = Coef_1*E + Coef_2*E_dot;

Coef_dot_1 = D_dot(1)/((1-D(1))*(1-D_dot(2)))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
Coef_dot_2 = 1/(1-D_dot(2))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));

d3_dot_inf = Coef_dot_1*E + Coef_dot_2*E_dot;

interval = (2*alpha+beta)/phi - d;
A1 = exp(-theta*interval)*cos(phi*interval) + theta/phi*exp(-theta*interval)*sin(phi*interval);
A2 = exp(-theta*interval)*sin(phi*interval)*1/phi;  
Bb = -(2*t*interval-1/theta*t*(1-k*t^2)*(exp(-theta*interval)*cos(phi*interval)-1)...
 -1/phi*t*(3-k*t^2)*exp(-theta*interval)*sin(phi*interval));
Bc = -(-t*(exp(-theta*interval)*cos(phi*interval)-1)- theta*t/phi*exp(-theta*interval)*sin(phi*interval));
B4= Bb*b(4)+Bc*c(4);

d4_inf =A1*d3_inf + A2*d3_dot_inf+B4;


D = simplify(D);
D_dot = simplify(D_dot);
A1 = simplify(A1);
A2 = simplify(A2);
E = simplify(E);
E_dot = simplify(E_dot);

%d3_inf = simplify(d3_inf)
d4_inf = simplify(d4_inf);

d4_inf_limit = limit(d4_inf,a,inf);

