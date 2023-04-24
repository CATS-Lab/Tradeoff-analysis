clear;

T =30;
tau = 0.9;
lambda = 0;
if lambda*tau>=1
    kappa = 2;
else
    %kappa = 0.9*(2-lambda*tau-2*sqrt(1-lambda*tau))/tau^2;
    kappa = 1
end

theta = (tau*kappa+lambda)/2;
varphi = sqrt((tau*kappa+lambda)^2-4*kappa)/2;
a_bar = 2;
a_underbar = 3; % maximum deceleration
v_bar = 20;



A = (theta-lambda-varphi)*(theta-varphi)*(exp((theta+varphi)*v_bar/a_underbar)-1);
B = (theta-lambda+varphi)*(theta+varphi)*(exp((theta-varphi)*v_bar/a_underbar)-1);
%t_11 = T-1/(2*varphi)*log(A/B)
t_11 =T-15;


TT=T+0.05;
ts = [0:0.01:TT];
xs = zeros(size(ts));
vs = zeros(size(ts));
ys = zeros(size(ts));
ws = zeros(size(ts));
as = zeros(size(ts));

for i = 2:length(ts)
    %update the first trajectory
    t = ts(i);
    delta = ts(i)-ts(i-1);
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
    
    as(i) = kappa*(xs(i-1)-ys(i-1)-tau*ws(i-1))+lambda*(vs(i-1)-ws(i-1));
    g = (as(i)-as(i-1))/delta;
    ws(i) = ws(i-1)+(as(i)+as(i-1))/2*delta;
    ys(i) = ys(i-1) + ws(i-1)*delta+as(i-1)*delta^2/2 ...
        +g*delta^3/6;
end

ys_ana = zeros(size(ts));
ws_ana = zeros(size(ts));


A = a_bar*(kappa*tau-(theta+varphi)*(kappa*tau^2+lambda*tau-1))/(2*varphi*kappa);
B = a_bar*(-kappa*tau+(theta-varphi)*(kappa*tau^2+lambda*tau-1))/(2*varphi*kappa);
C = a_bar/2;
D = -a_bar*tau;
E = a_bar*(kappa*tau^2+lambda*tau-1)/kappa;
for i = 1:length(ts)
    t = ts(i);
	if t<=v_bar/a_bar
        ys_ana(i)= A*exp((-theta+varphi)*t)+B*exp((-theta-varphi)*t)+C*t^2+D*t+E;
        ws_ana(i)=A*(-theta+varphi)*exp((-theta+varphi)*t)+B*(-theta-varphi)*exp((-theta-varphi)*t)+2*C*t+D;
    end
end

Hs = zeros (size(ts));

H_int = 0;
for i = 1:length(ts)
    t = ts(i);
    mu4_t = 1/(2*varphi)*(exp(-(theta+varphi)*(T-t))-exp(-(theta-varphi)*(T-t)));
    mu4_int = 1/(2*varphi)*((1-exp(-(theta+varphi)*(T-t)))/(theta+varphi)-(1-exp(-(theta-varphi)*(T-t)))/(theta-varphi));
    
    mu1_t = 1+kappa*mu4_int;
    Hs(i) = -mu1_t-lambda*mu4_t;
    if t>t_11 && t<=t_11+v_bar/a_underbar
         H_int = H_int + (Hs(i-1)+Hs(i))/2*(ts(i)-ts(i-1));
    end
end




figure(1);clf;hold all;
plot(ts,vs,'r');
plot(ts,ws,'b');
plot(ts,ws_ana,':k','linewidth',2);
figure(2);clf;hold all;
Diff = xs-ys;
L= length(ts)-1;
%L=10;
plot(ts(end-L:end),Diff(end-L:end),'r');
[Diff_min,i]=min(Diff);
plot(ts(i),Diff_min,'o');
Diff_min

Diff_min_ana = -A*exp((-theta+varphi)*v_bar/a_bar)-B*exp((-theta-varphi)*v_bar/a_bar)+tau*v_bar-a_bar*(kappa*tau^2+lambda*tau-1)/kappa;
plot(v_bar/a_bar,Diff_min_ana,'s');
title(num2str(Diff_min));

figure(3);clf;hold all;
plot(ts, Hs);
plot([t_11,t_11],[min(Hs),max(Hs)])

plot([t_11,t_11]+v_bar/a_underbar,[min(Hs),max(Hs)])

%H_int;
%H_max = max(Hs)
%[H_max,i] = max(Hs);
%ts(i)



