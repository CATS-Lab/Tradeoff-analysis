clear;

T =50;
tau = 0.5;
kappa = 0.9*4/tau^2;


a_bar = 2;
a_underbar = 3; % maximum deceleration
v_bar = 20;
t_11 =T-15;

omegas = [0:0.05:3];
ss = i*omegas;
TFs = (kappa)./(ss.^2+kappa+kappa*tau.*ss);
TF_max = max(abs(TFs));

TT=T+0.05;
delta = 0.01;


ts = [0: delta:TT];
xs = zeros(size(ts));
vs = zeros(size(ts));
ys = zeros(size(ts));
ws = zeros(size(ts));
as = zeros(size(ts));

for i = 2:length(ts)
    %update the first trajectory
    t = ts(i);
    if t<0
        continue
    end
   
    if t<=v_bar/a_bar
        vs(i) = min(vs(i-1)+a_bar*delta,v_bar);
    elseif t<=t_11
        vs(i) = vs(i-1);
    elseif t<=t_11+v_bar/a_underbar
        vs(i)=max(vs(i-1)-a_underbar*delta,0);
    else
        vs(i)=vs(i-1);
    end
    
    xs(i) = xs(i-1)+(vs(i)+vs(i-1))/2*delta;
    
    % carfollowing    
    as(i) = kappa*(xs(i-1)-ys(i-1)-tau*ws(i-1));
    g = (as(i)-as(i-1))/delta;
    ws(i) = ws(i-1)+(as(i)+as(i-1))/2*delta;
    ys(i) = ys(i-1) + ws(i-1)*delta+as(i-1)*delta^2/2 ...
        +g*delta^3/6;
end


figure(1);clf;hold all;
plot(ts,vs,'r');
plot(ts,ws,'b');

figure(4);clf;hold all;
Diff = xs-ys;
L= length(ts)-1;
%L=10;
plot(ts(end-L:end),Diff(end-L:end),'r');
[Diff_min,i]=min(Diff);
plot(ts(i),Diff_min,'o');
title(['min(Diff): ',num2str(Diff_min),'; TF_max: ',num2str(TF_max)])

Diff_min
TF_max
figure(3);clf;hold all;
plot(omegas ,abs(TFs),'b')
result = [tau,Diff_min,TF_max]

