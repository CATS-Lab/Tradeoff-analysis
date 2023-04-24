clear
k = [1:10000];

t = 4;
r = (-k.*t+sqrt(k.^2.*t.^2-4*k))/2;

figure(999);clf;hold on;
plot(k,r);