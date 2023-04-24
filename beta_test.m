clear
tau = 0.5;
kappa = 0.9*4/tau^2;

phi = sqrt(4*kappa-tau^2*kappa^2)/2;
theta = kappa*tau/2

deltas = [0:0.01:2*pi/phi];
%deltas = [3.1416    2.0944];

denom = exp(theta*deltas).*cos(phi*deltas)-1;
numer = exp(theta*deltas).*sin(phi*deltas);
betas = acot (denom./numer);

betas = betas+ pi.*((denom<0))+2*pi.*((denom>0).*(numer<0));

figure(1); clf; hold on;
plot(deltas*phi,betas);
    

figure (2);clf; hold on;
vec = [0.1:0.1:pi/2];
plot(vec,csc(vec));
%plot(deltas,(cos(phi*deltas)-exp(-theta*deltas))./sin(phi*deltas));