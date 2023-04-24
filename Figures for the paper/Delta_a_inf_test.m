clear
v = 2;
k = 2.2;
t = 0.21; 

phi = sqrt(4*k-t^2*k^2)/2;
theta = k*t/2;

%Delta_a_inf_0 = -((v*(exp((pi*k*t)/(-k*(k*t^2 - 4))^(1/2)) - 1)*(16*k^2*t*cos(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)) - 2*sin(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2))*(-k*(k*t^2 - 4))^(3/2) - 8*k^3*t^3*cos(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)) + k^4*t^5*cos(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)) + k*t^2*sin(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2))*(-k*(k*t^2 - 4))^(3/2)))/(k^2*(k*t^2 - 4)^2*(exp((k*t*(pi + atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)))/(-k*(k*t^2 - 4))^(1/2)) - 2*exp((k*t*(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)))/(-k*(k*t^2 - 4))^(1/2)) + exp((k*t*(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) - pi + 2*asin((k^(1/2)*t)/2)))/(-k*(k*t^2 - 4))^(1/2)))))
Delta_a_inf_0 = -((v*(exp((pi*k*t)/(-k*(k*t^2 - 4))^(1/2)) - 1)*(16*k^2*t*cos(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)) - 2*sin(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2))*(-k*(k*t^2 - 4))^(3/2) - 8*k^3*t^3*cos(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)) + k^4*t^5*cos(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)) + k*t^2*sin(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2))*(-k*(k*t^2 - 4))^(3/2)))/(k^2*(k*t^2 - 4)^2*(exp((k*t*(pi + atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)))/(-k*(k*t^2 - 4))^(1/2)) - 2*exp((k*t*(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) + 2*asin((k^(1/2)*t)/2)))/(-k*(k*t^2 - 4))^(1/2)) + exp((k*t*(atan((-k*(k*t^2 - 4))^(1/2)/(k*t)) - pi + 2*asin((k^(1/2)*t)/2)))/(-k*(k*t^2 - 4))^(1/2)))))


%Delta_a_inf = v*(exp(pi*theta/phi)-1) ...
        %/(k^(1/2)*exp((pi/2+atan(theta/phi))*theta/phi)*(exp(pi*theta/phi)+ exp(-pi*theta/phi)- 2))
        
%Delta_a_inf = v*(exp(pi*(4*k^-1*t^-2-1)^(-1/2))-1) ...
 %       /(k^(1/2)*exp((pi/2+atan(theta/phi))*theta/phi)*(exp(pi*theta/phi)+ exp(-pi*theta/phi)- 2))
 
 Delta_a_inf = v...
        /(k^(1/2)*exp(theta/phi*atan(theta/phi))*(exp(pi*theta/phi/2)- exp(-pi*theta/phi/2)))


Diff = Delta_a_inf_0-Delta_a_inf