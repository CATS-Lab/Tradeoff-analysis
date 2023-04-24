clear;
kappa = 0.2;
tau = 0.5*sqrt(4/kappa);




theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));


a_bar = 1;
a_underbar = a_bar; % maximum deceleration
v_bar = 10;

%delta_max = v_bar/a_bar; 
%ratio = min(2*pi/phi/(2*delta_max),1);
delta = min (v_bar/a_bar,pi/phi);
v_max = delta*a_underbar;


denom = exp(theta*delta).*cos(phi*delta)-1;
numer = exp(theta*delta).*sin(phi*delta);
beta = acot (denom./numer);
beta = beta+ pi.*(denom<0)+2*pi.*((denom>0).*(numer<0));

K=5;

T =(2*K-1)*pi/phi+(2*alpha+beta)/phi ;
%T =2;
%f_0 = mod(T-(alpha+pi/2)/phi,2*pi/phi)-2*pi/phi;
%K = (T-(alpha+pi/2)/phi-f_0)/(2*pi/phi);

t_dec_start = T-(2*alpha+beta)/phi-2*pi/phi*(K-[1:K]);
t_dec_end = t_dec_start+delta;

t_acc_start = t_dec_start - pi/phi;
t_acc_end = t_dec_end - pi/phi;
t_p = [t_acc_start;t_acc_end;t_dec_start;t_dec_end]; 
t_bar_p = t_p;

if size(t_bar_p,2)>0
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
end



TT=T+0.2;
interval = 0.001;


ts = [0:interval:TT];
xs = zeros(size(ts));
vs = zeros(size(ts));
ys = zeros(size(ts));
ws = zeros(size(ts));
as = zeros(size(ts));
ys_ana = zeros(size(ts));
ws_ana = zeros(size(ts));
diff = zeros(size(ts));
diff_ana = zeros(size(ts));
diff_dot_ana = zeros(size(ts));

j_old = 2;
y_0 = 0;
w_0 = 0;
t_0 = 0;
x_0 =0;
diff_0 = 0;
diff_dot_0 = 0;
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
    
    if ~is_break
        j =1;
    end
    
    if j==2
         a_temp = a_bar;
    elseif j==4
         a_temp = -a_bar;
    else
        a_temp =0;
    end
    
    
    vs(i) = vs(i-1)+a_temp*interval;
    %vs(i) = max(0,min(vs(i),v_bar));
    
    xs(i) = xs(i-1)+(vs(i)+vs(i-1))/2*interval;
    
    % carfollowing    
    as(i) = kappa*(xs(i-1)-ys(i-1)-tau*ws(i-1));
    g = (as(i)-as(i-1))/interval;
    ws(i) = ws(i-1)+(as(i)+as(i-1))/2*interval;
    ys(i) = ys(i-1) + ws(i-1)*interval+as(i-1)*interval^2/2 ...
        +g*interval^3/6;
    
    
    
    if j == 2
       b = 0.5* a_bar;
       c =0;       
    elseif j==3
       b = 0;
       c = delta* a_bar;       
    elseif j==4
       b = -0.5* a_bar;
       c = delta* a_bar;       
    else
        b = 0;
        c = 0;        
    end
    if j~=j_old
        x_0 = xs(i-1);
        y_0 = ys_ana(i-1)-x_0;
        w_0 = ws_ana(i-1);
        t_0 = ts(i-1);
        diff_0 = -y_0;
        diff_dot_0 = vs(i-1)-w_0  ;
       
    end
    
    t=ts(i)-t_0;
    
%     C = b;
%     D = c-2*tau*b;
%     A = y_0/2+C/kappa+tau*D/2;
%     B = 1/(2*phi)*(-2*theta*A+c-2*tau*b-w_0);
%     im = sqrt(-1);
%     ys_ana(i) = (A+im*B)*exp((-theta+im*phi)*t)+(A-im*B)*exp((-theta-im*phi)*t)+C*t^2+D*t+y_0-2*A;
    
    
   ys_ana(i) = (exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t))*y_0 + exp(-theta*t)*sin(phi*t)*1/phi*w_0...
        + (t^2-2*tau*t+2*tau/phi*exp(-theta*t)*sin(phi*t)+(1/theta*(exp(-theta*t)*cos(phi*t)-1)+1/phi*exp(-theta*t)*sin(phi*t))*(tau-2*theta*tau^2))*b...
        + (t+(exp(-theta*t)*cos(phi*t)-1)*tau-exp(-theta*t)*sin(phi*t)*(1-theta*tau)/phi)*c;
    
   D1 = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
    D2 = exp(-theta*t)*sin(phi*t)*1/phi;
    D3 = t^2-2*tau*t+1/theta*tau*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
        +1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
    D4 = t+(exp(-theta*t)*cos(phi*t)-1)*tau+ exp(-theta*t)*sin(phi*t)*(theta*tau-1)/phi;
    
    ys_ana(i) = D1*y_0 + D2*w_0 + D3*b+ D4*c;

    
    ys_ana(i) = ys_ana(i) + x_0;
    
    ws_ana(i) = -kappa/phi*exp(-theta*t)*sin(phi*t)*y_0+1/phi*(-theta*exp(-theta*t)*sin(phi*t)+phi*exp(-theta*t)*cos(phi*t))*w_0...
        +(2*t+2*tau*(exp(-theta*t)*cos(phi*t)-1)+1/phi*(2*theta*tau-2)*exp(-theta*t)*sin(phi*t))*b...
        +(1-exp(-theta*t)*cos(phi*t)-theta/phi*exp(-theta*t)*sin(phi*t))*c;
    
   
     A1 = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
     A2 = exp(-theta*t)*sin(phi*t)*1/phi;
     Bb = 2*tau*t-1/theta*tau*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
         -1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
     Bc = -tau*(exp(-theta*t)*cos(phi*t)-1)- theta*tau/phi*exp(-theta*t)*sin(phi*t);
    %C1 = D1;
    %C2 = D2;
    %C3 = -(D3-t^2);
    %C4 = -(D2+D4-t);
    
    
    diff_ana(i) = A1*diff_0+A2*diff_dot_0+Bb*b+Bc*c;
    diff (i) = xs(i)-ys(i);
    diff_dot(i) =(diff(i)-diff(i-1))/(ts(i)-ts(i-1));
    
     A1_dot = -kappa/phi*exp(-theta*t)*sin(phi*t);
     A2_dot  = -1/phi*theta*exp(-theta*t)*sin(phi*t)+exp(-theta*t)*cos(phi*t);
     Bb_dot = -2*tau*(exp(-theta*t)*cos(phi*t)-1)+1/phi*(2-2*theta*tau)*exp(-theta*t)*sin(phi*t);
     Bc_dot = tau*kappa/phi*exp(-theta*t)*sin(phi*t);
     diff_dot_ana (i)=  A1_dot*diff_0+A2_dot*diff_dot_0+Bb_dot*b+Bc_dot*c;
     j_old =j; 
end


%% Analytical solution
A = zeros (2,4);
B = zeros (1,4);
A_dot = zeros (2,4);
B_dot = zeros (1,4);

deltas = [delta,pi/phi-delta,delta,pi/phi-delta];
b = [0.5* a_bar,0,-0.5* a_bar,0];
c = [0,delta* a_bar,delta* a_bar,0];
 

for j = 1: 4
    t = deltas(j);
   
    A(1,j) = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
    A(2,j) = exp(-theta*t)*sin(phi*t)*1/phi;  
    Bb = 2*tau*t-1/theta*tau*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
     -1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
    Bc = -tau*(exp(-theta*t)*cos(phi*t)-1)- theta*tau/phi*exp(-theta*t)*sin(phi*t);
    B(j)= Bb*b(j)+Bc*c(j);
    
    A_dot(1,j) = -kappa/phi*exp(-theta*t)*sin(phi*t);
    A_dot(2,j)  = -1/phi*theta*exp(-theta*t)*sin(phi*t)+exp(-theta*t)*cos(phi*t);
    Bb_dot = -2*tau*(exp(-theta*t)*cos(phi*t)-1)+1/phi*(2-2*theta*tau)*exp(-theta*t)*sin(phi*t);
    Bc_dot = tau*kappa/phi*exp(-theta*t)*sin(phi*t);
    B_dot(j)= Bb_dot*b(j)+Bc_dot*c(j);
end
D = zeros(1,2);
D_dot = zeros(1,2);


E = B(3);
E_dot = B_dot(3);
for l3 = [1,2]
    D_temp3 = A(l3,3);
    D_dot_temp3 = A_dot(l3,3);    
    E_temp3 = A(l3,3);
    E_dot_temp3 = A_dot(l3,3);
    
    if l3 == 1
        E = E+E_temp3*B(2);
        E_dot = E_dot+E_dot_temp3*B(2);
    else
        E = E+E_temp3*B_dot(2);
        E_dot = E_dot+E_dot_temp3*B_dot(2);
    end
    for l2 = [1,2]
        if l3 == 1           
            D_temp2 = D_temp3*A(l2,2);
            D_dot_temp2 = D_dot_temp3*A(l2,2);            
            E_temp2 = E_temp3*A(l2,2);
            E_dot_temp2 = E_dot_temp3*A(l2,2);
        else
            D_temp2 = D_temp3 *A_dot(l2,2);
            D_dot_temp2 = D_dot_temp3 *A_dot(l2,2);
            E_temp2 = E_temp3*A_dot(l2,2);
            E_dot_temp2 = E_dot_temp3*A_dot(l2,2);
        end
        if l2 == 1
            E = E+E_temp2*B(1);
            E_dot = E_dot+E_dot_temp2*B(1);
        else
            E = E+E_temp2*B_dot(1);
            E_dot = E_dot+E_dot_temp2*B_dot(1);
        end
        for l1 = [1,2]
            if l2 == 1
                D_temp1 = D_temp2 *A(l1,1);
                D_dot_temp1 = D_dot_temp2 *A(l1,1);
                E_temp1 = E_temp2 *A(l1,1);
                E_dot_temp1 = E_dot_temp2 *A(l1,1);
            else
                D_temp1 = D_temp2 *A_dot(l1,1);
                D_dot_temp1 = D_dot_temp2 *A_dot(l1,1);
                E_temp1 = E_temp2 *A_dot(l1,1);
                E_dot_temp1 = E_dot_temp2 *A_dot(l1,1);
            end
            
            if l1 ==1
                D = D+ D_temp1 *[A(1,4),A(2,4)];
                D_dot = D_dot+ D_dot_temp1 *[A(1,4),A(2,4)];
                E = E+E_temp1*B(4);
                E_dot = E_dot+E_dot_temp1*B(4);
            else
                D = D+ D_temp1 *[A_dot(1,4),A_dot(2,4)];
                D_dot = D_dot+ D_dot_temp1 *[A_dot(1,4),A_dot(2,4)];
                E = E+E_temp1*B_dot(4);
                E_dot = E_dot+E_dot_temp1*B_dot(4);
            end
        end
    end
end

ks = [1:K];
t3 = -pi/phi + delta + ks*2*pi/phi;
Delta3 = E;
Delta_dot_3 = E_dot;
for k = 2:K
    Delta3(k) = D(1)*Delta3(k-1)+D(2)*Delta_dot_3(k-1)+E;
    Delta_dot_3(k) = D_dot(1)*Delta3(k-1)+D_dot(2)*Delta_dot_3(k-1)+E_dot;
end



%Coef_1 =  1/(1-D(1))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
%Coef_2 = D(2)/((1-D(1))*(1-D_dot(2)))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
%Delta3_inf = Coef_1*E + Coef_2*E_dot;
D(1) = exp(-2*pi*kappa*tau/(sqrt(4*kappa-tau^2*kappa^2)));
Delta3_inf  = 1/(1-D(1))*E;

% Coef_dot_1 = D_dot(1)/((1-D(1))*(1-D_dot(2)))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
% Coef_dot_2 = 1/(1-D_dot(2))*1/(1-D_dot(1)*D(2)/((1-D(1))*(1-D_dot(2))));
% 
% Delta3_dot_inf = Coef_dot_1*E + Coef_dot_2*E_dot;
D_dot(2)=  exp(-2*pi*kappa*tau/(sqrt(4*kappa-tau^2*kappa^2)));
Delta3_dot_inf  = 1/(1-D_dot(2))*E_dot;

t = (2*alpha+beta)/phi - delta;
A1 = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
A2 = exp(-theta*t)*sin(phi*t)*1/phi;  
Bb = 2*tau*t-1/theta*tau*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
 -1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
Bc = -tau*(exp(-theta*t)*cos(phi*t)-1)- theta*tau/phi*exp(-theta*t)*sin(phi*t);
B4= Bb*b(4)+Bc*c(4);

Delta4_inf =A1*Delta3_inf + A2*Delta3_dot_inf+B4;


gap =min(diff)-Delta4_inf;

%% Plot

figure(1);
subplot(2,1,1);cla;hold on;

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
if size(t_bar_p,2)>0
    sep = t_p(2,1)-t_p(1,1);
    plot([t_p(1,1),t_p(2,1)],[0,sep*a_bar],'--r');
    plot([t_p(2,1),t_p(3,1)],[sep*a_bar,sep*a_bar],'--r')
    plot([t_p(3,1),t_p(4,1)],[sep*a_bar,0],'--r')
    xLimits = get(gca,'XLim');
end

subplot(2,1,2);cla;hold on;
t = [0:0.1:T];
H_var = cos(phi*(T-t)-(alpha+pi)); %H /exp(-theta*(T-t))/sqrt(kappa)/phi
plot(t,H_var);

plot([0,0],[-1,1],':k');
plot([0,T],[0,0],'--k');
if size(t_bar_p,2)>0
    plot([t_bar_p(2,1),t_bar_p(2,1)],[-1,1],'g');
    plot([t_bar_p(4,1),t_bar_p(4,1)],[-1,1],'r');
    %plot([t_1,t_1],[-1,1],'k');
    %plot([t_2,t_2],[-1,1],'r');
    
    xlim(xLimits);
end



figure(2);
subplot(2,1,1);cla;hold all;
plot(ts,vs,'b');
plot(ts,ws,'r','linewidth',2);
plot(ts,ws_ana,'y--','linewidth',2);
plot(T,0,'o');

subplot(2,1,2);cla;hold all;
plot(ts,xs,'b');
plot(ts,ys,'r','linewidth',2);
plot(ts,ys_ana,'y--','linewidth',2);

figure(3);
subplot(2,1,1);cla;hold all;
%plot(ts,vs,'b');

for t = t_bar_p(4,:)
    plot(t*[1,1],[min(diff),max(diff)],'k--');
end
plot(ts,diff,'r','linewidth',1);
plot(ts,diff_ana,'g--','linewidth',1);
[diff_min,i]=min(diff_ana);
plot(T,diff_min,'o');
title(['min(Diff): ',num2str(diff_min)])

plot(t3, Delta3,'*');
plot([min(ts),max(ts)],Delta3_inf*[1,1],'k--');
plot([0,T],[Delta4_inf,Delta4_inf],'b:*');


subplot(2,1,2);cla;hold all;
for t = t_bar_p(4,:)
    plot(t*[1,1],[min(diff_dot),max(diff_dot)],'k--');
end
plot(ts,diff_dot,'r','linewidth',1);
plot(ts,diff_dot_ana,'g--','linewidth',1);
plot(T,0,'o');
[diff_dot_min,i]=min(diff_dot_ana);
plot(t3, Delta_dot_3,'*');
plot([min(ts),max(ts)],Delta3_dot_inf*[1,1],'k--');
title(['min(Diff_dot): ',num2str(diff_dot_min)])

diff_min;
diff_dot_min;
E
E_dot
Delta4_inf



% omegas = [0:0.01:3];
% i = sqrt(-1);
% ss = i*omegas;
% TFs = (kappa)./(ss.^2+kappa+kappa*tau.*ss);
% TF_max = max(abs(TFs));
% 
% 
% figure(4);clf;hold all;
% plot(omegas ,abs(TFs),'b')
% result = [tau,Diff_min,TF_max];
% title(['TF_max: ',num2str(TF_max)])

