clear;

T =12;
kappa = 1;
tau = sqrt(0.5/kappa);


theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));


a_bar = 20;
a_underbar = a_bar; % maximum deceleration
v_bar = 20;

delta_max = v_bar/a_bar; 
ratio = min(2*pi/phi/(2*delta_max),1);
delta = delta_max*ratio;
v_max = delta*a_underbar;


denom = exp(theta*delta).*cos(phi*delta)-1;
numer = exp(theta*delta).*sin(phi*delta);
beta = acot (denom./numer);
beta = beta+ pi.*(denom<0)+2*pi.*((denom>0).*(numer<0));


t_0e = T-(alpha+pi)/phi;
phase_0 = 2*pi/phi- mod(t_0e,2*pi/phi);

t_dec_1e = T-(2*alpha+beta)/phi;
t_acc_1e = T-(2*alpha+beta)/phi+pi/phi;
phase_acc_1 = t_acc_1e-t_0e;
v_max_gap = pi/phi - delta ;

gamma = mod(T-(alpha+pi)/phi,2*pi/phi) - 2*pi/phi;

k = (T-(alpha+pi)/phi - gamma)/(2*pi/phi);

%t_0 = mod(T-alpha/phi-pi/2/phi,2*pi/phi)
%t_dec_11 = mod(t_dec_1e,2*pi/phi);
t_dec_111 = T-(2*alpha+beta)/phi-2*k*pi/phi;
%t_dec_21 = mod(t_dec_1e+delta,2*pi/phi);
t_dec_211 = t_dec_111+delta;
%t_acc_11 = mod(t_dec_1e-pi/phi,2*pi/phi);
t_acc_111 = t_dec_111 +pi/phi;
%t_acc_21 = mod(t_dec_1e+delta-pi/phi,2*pi/phi);
t_acc_211 = t_dec_211 +pi/phi;
t_crit = [t_dec_1e,t_dec_1e+delta, T];
v_crit = [v_max,0,0];

while 1
    t_next = t_crit(1:2)-pi/phi;
    
    if t_next(2)>=0
        if t_next(2)<t_crit(1)
            v_next = v_crit(1);
            t_crit = [t_next(2),t_crit];
            v_crit = [v_next,v_crit];
        end
    else
        break;
    end
    
    if t_next(1)>=0
        if t_next(1)<t_crit(1)
            v_next = v_max-v_crit(1);
            t_crit = [t_next(1),t_crit];
            v_crit = [v_next,v_crit];
        end
    else
        break;
    end
    
end
% if v_crit(1)==0
%     t_crit = [0,t_crit];
%     v_crit = [0,v_crit];
% end

if phase_0>=0.25*2*pi/phi && phase_0<=phase_acc_1
   if 0<t_crit(1)
        t_crit = [0,t_crit];
        v_crit = [0,v_crit];
   end
% elseif phase_0 > phase_acc_1 &&  phase_0 < phase_acc_1 + v_max_gap
%     if delta> t_crit(1)
%         t_crit(1)=delta;
%     else
%         t_crit = [delta,t_crit];
%         v_crit = [v_max,v_crit];
%     end
%     
%     t_crit = [0,t_crit];
%     v_crit = [0,v_crit];   
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
    
    if delta_1> delta_max
         delta_1 = delta_max;
         if delta_1<t_crit(1)
             t_crit = [delta_1,t_crit];
             v_crit =  [v_bar,v_crit];
         elseif delta_1>t_crit(1)
             t_crit(1)=delta_1;
         end
         if 0<t_crit(1)
             t_crit = [delta_1,t_crit];
             v_crit =  [v_bar,v_crit];
         end
         t_crit = [0,t_crit];
         v_crit =  [0,v_crit];
    else
        i_next = find(v_max-v_crit,1);
        t_crit = t_crit(i_next:end);
        v_crit = v_crit(i_next:end);
        if 0<t_crit(1)
            t_crit = [0,delta_1,2*delta_1, t_crit];
            v_crit =  [0,delta_1*a_bar,0,v_crit];
        end
    end
    
end
    
figure(1);clf;
subplot(2,1,1);hold on;
plot(t_crit,v_crit);
%plot([t_1,t_1],[0,v_max],'k');
%plot([t_2,t_2],[0,v_max],'r');

plot([t_dec_211,t_dec_211],[0,v_max],'b','linewidth',1);
plot([t_acc_211,t_acc_211],[0,v_max],'r','linewidth',1);
plot([t_dec_111,t_dec_111],[0,v_max],'b+','linewidth',2);
plot([t_acc_111,t_acc_111],[0,v_max],'r+','linewidth',2);


subplot(2,1,2);hold on;
t = [0:0.1:T];
H_var = cos(phi*(T-t)-(alpha+pi)); %H /exp(-theta*(T-t))/sqrt(kappa)/phi
plot(t,H_var);

t_temp = mod((T- 2*alpha/phi),pi/phi);
%plot([t_temp,t_temp],[-1,1]);

plot([t_crit(2),t_crit(2)],[-1,1],':k');
plot([t_crit(3),t_crit(3)],[-1,1],':r');
%plot([t_1,t_1],[-1,1],'k');
%plot([t_2,t_2],[-1,1],'r');
plot([0,T],[0,0],'--k');







TT=T+0.05;
interval = 0.0001;


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
    
    for j = 2:length(t_crit)
        if t<=t_crit(j)
            break;
        end
    end
    a_temp = (v_crit(j)-v_crit(j-1))/(t_crit(j)- t_crit(j-1));
        
    
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

Diff_min


omegas = [0:0.1:3];
i = sqrt(-1);
ss = i*omegas;
TFs = (kappa)./(ss.^2+kappa+kappa*tau.*ss);
TF_max = max(abs(TFs));


figure(4);clf;hold all;
plot(omegas ,abs(TFs),'b')
result = [tau,Diff_min,TF_max]
title(['TF_max: ',num2str(TF_max)])

