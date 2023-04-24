kappa = 0.01;
tau = 1.5; %sec



delta = 0.1; % time increases interval sec
T = 1000;
ts = [0:delta:T];
N = length (ts);
K = 0; %Zero padding rate 

%% analytical solution
omegas = [0:N*(1+K)-1]*2*pi/N; %digital angular frequency 
Omegas = omegas./delta; %analog angular frequency sec^-1

TF = kappa./(-Omegas.^2+kappa+1i*kappa*tau*Omegas);


%% numerical solution

vs1 =  zeros(1,N);
vs1(1)=10;
xs1 = zeros(size(vs1));
for n = 2:N
    xs1(n) = xs1(n-1) + vs1(n-1)*delta; 
end

as2 = zeros(size(vs1));
vs2 = zeros(size(vs1));
xs2 = zeros(size(xs1));
for n = 1:N-1
    as2(n)= kappa*(xs1(n)-xs2(n)-tau*vs2(n));
    vs2(n+1)= vs2(n)+ as2(n)*delta;
    xs2(n+1) = xs2(n) + (vs2(n)+vs2(n+1))/2*delta;
end

Vs1 = fft([vs1,zeros(1,K*N)]);
Vs2 = fft([vs2,zeros(1,K*N)]);
TF_num = Vs2./Vs1;

%% plots

figure(1);clf;hold all;
plot(ts,xs1,'b');
plot(ts,xs2,'--r');

    


figure(2);clf;hold all;
plot(Omegas,abs(TF),'b');
plot(Omegas/(K+1),abs(TF_num),'--r');
 xlim([0,1])   

