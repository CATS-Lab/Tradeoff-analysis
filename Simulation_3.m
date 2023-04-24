clear;

T =50;
tau = 0.5;
kappa = 0.9*4/tau^2;

theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));


a_bar = 100;
a_underbar = 1; % maximum deceleration
v_bar = 20;
delta_acc_max = v_bar/a_bar; 
delta_dec_max = v_bar/a_underbar;
ratio = 2*pi/phi/(delta_acc_max+delta_dec_max);
delta(1) = delta_acc_max*ratio;
delta(2) = delta_dec_max*ratio;


denom = exp(theta*delta).*cos(phi*delta)-1;
numer = exp(theta*delta).*sin(phi*delta);
beta = acot (denom./numer);
beta = beta+ pi.*(denom<0)+2*pi.*((denom>0).*(numer<0));

phase_e = T-(alpha+pi)/phi;

t1e(1) = (T-(2*alpha+beta(2))/phi);
t1e(2) = t1e(1)+delta(2);

(t1e-phase_e)/ (pi/phi) 

return





TT=T+0.05;
interval = 0.01;


ts = [0: interval:TT];
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
        vs(i) = min(vs(i-1)+a_bar*interval,v_bar);
    elseif t<=t_11
        vs(i) = vs(i-1);
    elseif t<=t_11+v_bar/a_underbar
        vs(i)=max(vs(i-1)-a_underbar*interval,0);
    else
        vs(i)=vs(i-1);
    end
    
    xs(i) = xs(i-1)+(vs(i)+vs(i-1))/2*interval;
    
    % carfollowing    
    as(i) = kappa*(xs(i-1)-ys(i-1)-tau*ws(i-1));
    g = (as(i)-as(i-1))/interval;
    ws(i) = ws(i-1)+(as(i)+as(i-1))/2*interval;
    ys(i) = ys(i-1) + ws(i-1)*interval+as(i-1)*interval^2/2 ...
        +g*interval^3/6;
end


figure(1);clf;hold all;
plot(ts,vs,'r');
plot(ts,ws,'b');


omegas = [0:0.05:3];
ss = i*omegas;
TFs = (kappa)./(ss.^2+kappa+kappa*tau.*ss);
TF_max = max(abs(TFs));

figure(2);clf;hold all;
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

