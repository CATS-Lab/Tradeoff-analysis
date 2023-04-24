clear;

T =16.5;
kappa = 1;
tau = sqrt(0.5/kappa);


theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));


a_bar = 20;
a_underbar = a_bar; % maximum deceleration
v_bar = 20;

%delta_max = v_bar/a_bar; 
%ratio = min(2*pi/phi/(2*delta_max),1);
delta = min (v_bar/a_bar,pi/phi);
v_max = delta*a_underbar;


denom = exp(theta*delta).*cos(phi*delta)-1;
numer = exp(theta*delta).*sin(phi*delta);
beta = acot (denom./numer);
beta = beta+ pi.*(denom<0)+2*pi.*((denom>0).*(numer<0));


f_0 = mod(T-(alpha+pi/2)/phi,2*pi/phi)-2*pi/phi;
K = (T-(alpha+pi/2)/phi-f_0)/(2*pi/phi);

t_dec_start = T-(2*alpha+beta)/phi-2*pi/phi*(K-[1:K]);
t_dec_end = t_dec_start+delta;

t_acc_start = t_dec_start - pi/phi;
t_acc_end = t_dec_end - pi/phi;
t_p = [t_acc_start;t_acc_end;t_dec_start;t_dec_end]; 
t_bar_p = t_p;

if t_p(3,1)-delta>=0
    t_bar_p(1,1)= max(0,t_p(1,1));
    t_bar_p(2,1)= t_bar_p(1,1)+delta;
else
    t_temp = mod(T- 2*alpha/phi,2*pi/phi);
    old_lhs = -t_temp;
    inc = 0.001;
    for delta1 = [inc:inc:2*pi/phi]
        denom1 = exp(theta*delta1).*cos(phi*delta1)-1;
        numer1 = exp(theta*delta1).*sin(phi*delta1);
        beta1 = acot (denom1/numer1);
        beta1 = beta1+ pi.*(denom1<0)+2*pi.*((denom1>0).*(numer1<0));
        new_lhs = beta1/phi+delta1-t_temp;
        if new_lhs>=0
            delta_1 = delta1-new_lhs/(new_lhs-old_lhs)*inc;
            break
        end        
        old_lhs = new_lhs;
    end 
    t_bar_p(1,1)= 0;
    t_bar_p(2,1)= delta_1;
    t_bar_p(3,1)= delta_1;
    t_bar_p(4,1)= 2*delta_1;
    
end

    
figure(1);clf;
subplot(2,1,1);hold on;

for i =1:size(t_bar_p,2)
    sep = t_bar_p(2,i)-t_bar_p(1,i);
    if i ==1
        plot([0,t_bar_p(1,i)],[0,0],'b');
    else
        plot([t_bar_p(4,i-1),t_bar_p(1,i)],[0,0],'b');
    end
    plot([t_bar_p(1,i),t_bar_p(2,i)],[0,sep*a_bar],'b');
    plot([t_bar_p(2,i),t_bar_p(3,i)],[sep*a_bar,sep*a_bar],'b')
    plot([t_bar_p(3,i),t_bar_p(4,i)],[sep*a_bar,0],'b')
    
    if i ==size(t_p,2)
        plot([t_bar_p(4,i),T],[0,0],'b');
    end
end

sep = t_p(2,1)-t_p(1,1);
plot([t_p(1,1),t_p(2,1)],[0,sep*a_bar],'--r');
plot([t_p(2,1),t_p(3,1)],[sep*a_bar,sep*a_bar],'--r')
plot([t_p(3,1),t_p(4,1)],[sep*a_bar,0],'--r')
xLimits = get(gca,'XLim');

subplot(2,1,2);hold on;
t = [0:0.1:T];
H_var = cos(phi*(T-t)-(alpha+pi)); %H /exp(-theta*(T-t))/sqrt(kappa)/phi
plot(t,H_var);

plot([0,0],[-1,1],':k');
plot([t_bar_p(2,1),t_bar_p(2,1)],[-1,1],'g');
plot([t_bar_p(4,1),t_bar_p(4,1)],[-1,1],'r');
%plot([t_1,t_1],[-1,1],'k');
%plot([t_2,t_2],[-1,1],'r');
plot([0,T],[0,0],'--k');
xlim(xLimits);






TT=T+0.2;
interval = 0.00001;


ts = [0:interval:TT];
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
    
    is_break =0;
    for k = 1:size(t_bar_p,2)
        for j = 1:4
            if t<=t_bar_p(j,k)
                is_break =1;
                break;
            end
        end
        if is_break
            break;
        end
    end
    
    if j==2
         a_temp = a_bar;
    elseif j==4
         a_temp = -a_bar;
    else
        a_temp =0;
    end
    
    
    vs(i) = vs(i-1)+a_temp*interval;
    vs(i) = max(0,min(vs(i),v_bar));
    
    xs(i) = xs(i-1)+(vs(i)+vs(i-1))/2*interval;
    
    % carfollowing    
    as(i) = kappa*(xs(i-1)-ys(i-1)-tau*ws(i-1));
    g = (as(i)-as(i-1))/interval;
    ws(i) = ws(i-1)+(as(i)+as(i-1))/2*interval;
    ys(i) = ys(i-1) + ws(i-1)*interval+as(i-1)*interval^2/2 ...
        +g*interval^3/6;
end


figure(2);clf;hold all;
plot(ts,vs,'r');
plot(ts,ws,'b');




figure(3);clf;hold all;
Diff = xs-ys;
L= length(ts)-1;
%L=10;
plot(ts(end-L:end),Diff(end-L:end),'r');
[Diff_min,i]=min(Diff);
plot(T,Diff_min,'o');
title(['min(Diff): ',num2str(Diff_min)])

ts(i)
Diff_min


omegas = [0:0.1:3];
i = sqrt(-1);
ss = i*omegas;
TFs = (kappa)./(ss.^2+kappa+kappa*tau.*ss);
TF_max = max(abs(TFs));


figure(4);clf;hold all;
plot(omegas ,abs(TFs),'b')
result = [tau,Diff_min,TF_max];
title(['TF_max: ',num2str(TF_max)])

