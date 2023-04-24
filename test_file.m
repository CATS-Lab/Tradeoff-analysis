% A =phi* (kappa*tau^2-1);
% B = theta*(kappa*tau^2-3);
clear
% C = kappa;
% A^2+B^2-C
% 
% D = tau*phi;
% E = (1-tau*theta);
% 
% D^2+E^2-1
% 
% ts = [0:0.05:100];
% ys = 3*cos(ts) +5*cos(ts+2);
% figure(10);cla; hold all;
% plot(ts,ys);
% 
% 
% H = -delta-0.5/theta*(kappa*tau^2-1);
% I = -delta*theta/phi -0.5/phi *(kappa*tau^2-3);
% J = delta^2 *kappa/phi^2 +1/(2*theta*tau*phi^2)+delta/(theta*phi^2)*(2*theta^2-kappa);
% H^2+I^2-J

% w = [-9:0.1:9];
% y = acot(w);
% figure(123);cla;
% plot(w,y)
int =0.001;
x=[int:int:1];
denom = cos(x)-1;
numer = sin(x);
y = atan(numer./denom)+;
