clear;
kappa = 0.2;
tau = 0.2*sqrt(0.5/kappa);

theta = (tau*kappa)/2;
phi = sqrt(4*kappa-(tau*kappa)^2)/2;
alpha = asin(theta/sqrt(kappa));


for iter = 1:3
    if iter ==1
        T=45;
        a_bar = 5;
    elseif iter ==2
        T=45;
        a_bar = 2.8;
    elseif iter ==3
        T=37;
        a_bar = 5;
        
    end
    interval = 0.1;
    ts = [0: interval:T];
    nmu1 = sqrt(kappa)/phi*exp(-theta*(T-ts)).*sin(phi*(T-ts)-alpha-pi/2);
    
    
    
    v_bar = 20;

    delta = min (v_bar/a_bar,pi/phi);
    v_max = delta*a_bar;

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

    t1 =0;
    t2 =1;

    mu2 = zeros(size(ts));
    for i = length(ts)-1:-1:1
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

        if ((j==2 || j==4) && is_break)|| (k ==1 && j==1)
            %a_temp = a_bar; 

            t2 = t;
            if k ==1 && j<=2
                t1 = t_bar_p(2,k);
            else
                t1 = max(0,t_bar_p(j-1,k));
            end
            t2 = t; 
            mu2(i) = (kappa^(1/2)*((exp(-(kappa*tau*(T - t1))/2)*(kappa^(3/2)*tau^2*sin(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t1))/2)...
                - 2*sin(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t1))/2)*(- kappa^2*tau^2 + 4*kappa)^(1/2)*(1 - (kappa*tau^2)/4)^(1/2)...
                + kappa^(1/2)*tau*cos(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t1))/2)*(- kappa^2*tau^2 + 4*kappa)^(1/2)))/(2*kappa) ...
                - (exp(-(kappa*tau*(T - t2))/2)*(kappa^(3/2)*tau^2*sin(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t2))/2) ...
                - 2*sin(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t2))/2)*(- kappa^2*tau^2 + 4*kappa)^(1/2)*(1 - (kappa*tau^2)/4)^(1/2) ...
                + kappa^(1/2)*tau*cos(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t2))/2)*(- kappa^2*tau^2 + 4*kappa)^(1/2)))/(2*kappa)...
                + tau*exp(-(kappa*tau*(T - t1))/2)*cos(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t1))/2)*(1 - (kappa*tau^2)/4)^(1/2)...
                - tau*exp(-(kappa*tau*(T - t2))/2)*cos(((- kappa^2*tau^2 + 4*kappa)^(1/2)*(T - t2))/2)*(1 - (kappa*tau^2)/4)^(1/2)))/(4*kappa - kappa^2*tau^2)^(1/2);


        else

            mu2(i)=0;

        end
    end
    ,



    figure(9+iter);clf;hold all;
    
     set(gca,'TickLabelInterpreter','latex');
    if iter ==1
        set(gcf,'Position', [300,880-250*iter,600,150]);
        set(gca,'InnerPosition',[0.03 0.21 0.95 0.7]);
    else
        set(gcf,'Position', [300,880-250*iter,600,130]);
        set(gca,'InnerPosition',[0.03 0.22 0.95 0.63]);
    end
    
    %box on;
    plot(ts,nmu1,':m','linewidth',2);
    plot(ts,mu2*.6,'r','linewidth',2);
    


    label = '--b';
    lw = 2;

    for i =1:size(t_bar_p,2)
        sep = (t_bar_p(2,i)-t_bar_p(1,i))/15;
        if i ==1
            plot([0,t_bar_p(1,i)],[0,0],label,'linewidth',lw);
        else
            plot([t_bar_p(4,i-1),t_bar_p(1,i)],[0,0],label,'linewidth',lw);
        end
        plot([t_bar_p(1,i),t_bar_p(2,i)],[0,sep*a_bar],label,'linewidth',lw);
        plot([t_bar_p(2,i),t_bar_p(3,i)],[sep*a_bar,sep*a_bar],label,'linewidth',lw)
        plot([t_bar_p(3,i),t_bar_p(4,i)],[sep*a_bar,0],label,'linewidth',lw)

        if i ==size(t_p,2)
            plot([t_bar_p(4,i),T],[0,0],'b');
        end
    end

    plot([ts(1),ts(end)],[0,0],':k');
    sep1 = (t_bar_p(2,1)-t_bar_p(1,1))/15;
    sepe = (t_bar_p(2,end)-t_bar_p(1,end))/15;
    
    if iter <=2
        plot([ts(1),t_bar_p(2,1)],[sep1*a_bar,sep1*a_bar],':k')
    else
        plot([ts(1),t_bar_p(2,2)],[sepe*a_bar,sepe*a_bar],':k')
        
    end
    
    yl = ylim;
    if iter ==1
        ylim([-1,yl(2)*1.2]);
    else
        ylim([-1,yl(2)]);
    end
    
    yticks([-1,0,sepe*a_bar]);
    yticklabels({'-1','0','$\bar{v}$'});

   
    if iter ==1
        title('(a)')
    elseif iter ==2
        title('(b)')
    elseif iter ==3
        title('(c)')
    end
    
    set(get(gca,'title'),'Position',[2, 1.5 .5])

 
   



    plot([t_bar_p(1,1),t_bar_p(1,1)],[yl(1),0],':k');
    plot([t_bar_p(2,1),t_bar_p(2,1)],[yl(1),sep1*a_bar],':k');
    plot([t_bar_p(3,1),t_bar_p(3,1)],[yl(1),sep1*a_bar],':k');
    plot([t_bar_p(4,1),t_bar_p(4,1)],[yl(1),0],':k');

    plot([t_bar_p(1,end),t_bar_p(1,end)],[yl(1),0],':k');
    plot([t_bar_p(2,end),t_bar_p(2,end)],[yl(1),sepe*a_bar],':k');
    plot([t_bar_p(3,end),t_bar_p(3,end)],[yl(1),sepe*a_bar],':k');
    plot([t_bar_p(4,end),t_bar_p(4,end)],[yl(1),0],':k');
    
  

    if iter ==1
    xticks([0,t_bar_p(1,1),t_bar_p(2,1),t_bar_p(3,1),t_bar_p(4,1),...
        t_bar_p(1,end),t_bar_p(2,end),t_bar_p(3,end),t_bar_p(4,end),ts(end)]);
    xticklabels({'0','$t_1^{\mbox{acc}-}$','$t_1^{\mbox{acc}+}$','$t_1^{\mbox{dec}-}$','$t_1^{\mbox{dec}+}$',...
       '$t_K^{\mbox{acc}-}$','$t_K^{\mbox{acc}+}$','$t_K^{\mbox{dec}-}$','$t_K^{\mbox{dec}+}$','$T$'});
    elseif iter ==2
        xticks([0,t_bar_p(1,1),t_bar_p(2,1),t_bar_p(4,1),...
        t_bar_p(1,end),t_bar_p(2,end),ts(end)]);
        xticklabels({'0','$t_1^{\mbox{acc}-}$','$t_1^{\mbox{acc}+}/t_1^{\mbox{dec}-}$','$t_1^{\mbox{dec}+}/t_2^{\mbox{acc}-}$',...
       '$t_K^{\mbox{acc}-}/t_{K-1}^{\mbox{dec}+}$','$t_K^{\mbox{acc}+}/t_K^{\mbox{dec}-}$','$T$'});
    elseif iter ==3
         xticks([0,t_bar_p(3,1),t_bar_p(4,1),...
       t_bar_p(1,end),t_bar_p(2,end),t_bar_p(3,end),t_bar_p(4,end),ts(end)]);
        xticklabels({'$0$','$t_1^{\mbox{acc}+}/t_1^{\mbox{dec}-}$','$t_1^{\mbox{dec}+}$',...
       '$t_K^{\mbox{acc}-}$','$t_K^{\mbox{acc}+}$','$t_K^{\mbox{dec}-}$','$t_K^{\mbox{dec}+}$','$T$'});

    end
    
    if iter<=2
        x_temps = [(t_bar_p(1,1)+t_bar_p(2,1))/2,(t_bar_p(3,1)+t_bar_p(4,1))/2];
        for x_temp = x_temps
            text(x_temp,-0.8,'$\bar{\delta}$','interpreter','latex');
            %set(tx,'interpreter','latex');
        end
    else
        x_temps = [(t_bar_p(1,1)+t_bar_p(2,1))/2,(t_bar_p(3,1)+t_bar_p(4,1))/2];
        for x_temp = x_temps
            text(x_temp,-0.8,'$\delta_1$','interpreter','latex');
            %set(tx,'interpreter','latex');
        end
    end
    x_temps = [(t_bar_p(1,end)+t_bar_p(2,end))/2,(t_bar_p(3,end)+t_bar_p(4,end))/2];
    for x_temp = x_temps
        text(x_temp,-0.8,'$\bar{\delta}$','interpreter','latex');
        %set(tx,'interpreter','latex');
    end
   
    xlabel('$t$','interpreter','latex')

    if iter ==1
        lh=legend({'$-\mu_1(t)$','$\mu_2(t)$','$v(t)$'},'interpreter','latex','Position',[0.1 0.85 0.8 0.1]);
        set(lh,'Orientation','horizontal')
        set(lh,'Box','off')
    end
    saveas(gcf,['../../Figures/illurstration_mu_v_',num2str(iter)],'png');
end
