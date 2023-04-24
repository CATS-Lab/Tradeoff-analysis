clear;

kappas = [0.5,0.5,1,1];
taus = [0.5,1,0.5,1];
v_bar = 1;

subfig_ids = {'(a)','(b)','(c)','(d)'};
for iter = 1:length(kappas)
    kappa = kappas(iter);
    tau = taus(iter);
    rs = [0:0.1:5];
    Deltas=zeros(size(rs));
    
    for i_vec = 1:length(rs)
       
        a_bar = rs(i_vec)*v_bar;
        
        
        if tau^2*kappa>=4
            continue
        end


        theta = (tau*kappa)/2;
        phi = sqrt(4*kappa-(tau*kappa)^2)/2;
        alpha = asin(theta/sqrt(kappa));
        delta = min (v_bar/a_bar,pi/phi);
        v_max = delta*a_bar;

        denom = exp(theta*delta).*cos(phi*delta)-1;
        numer = exp(theta*delta).*sin(phi*delta);
        beta = acot (denom./numer);
        beta = beta+ pi.*(denom<0)+2*pi.*((denom>0).*(numer<0));


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
            Bb = -2*tau*t-2/kappa*(kappa*tau^2-1)*(exp(-theta*t)*cos(phi*t)-1)...
             -1/phi*tau*(kappa*tau^2-3)*exp(-theta*t)*sin(phi*t);
            Bc = tau*(exp(-theta*t)*cos(phi*t)-1)+ theta*tau/phi*exp(-theta*t)*sin(phi*t);
            B(j)= Bb*b(j)+Bc*c(j);

            A_dot(1,j) = -kappa/phi*exp(-theta*t)*sin(phi*t);
            A_dot(2,j)  = -1/phi*theta*exp(-theta*t)*sin(phi*t)+exp(-theta*t)*cos(phi*t);
            Bb_dot = 2*tau*(exp(-theta*t)*cos(phi*t)-1)-1/phi*(2-2*theta*tau)*exp(-theta*t)*sin(phi*t);
            Bc_dot = -tau*kappa/phi*exp(-theta*t)*sin(phi*t);
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




        D(1) = exp(-2*pi*kappa*tau/(sqrt(4*kappa-tau^2*kappa^2)));
        Delta3_inf  = 1/(1-D(1))*E;

        D_dot(2)=  exp(-2*pi*kappa*tau/(sqrt(4*kappa-tau^2*kappa^2)));
        Delta3_dot_inf  = 1/(1-D_dot(2))*E_dot;

        t = (2*alpha+beta)/phi - delta;
        A1 = exp(-theta*t)*cos(phi*t) + theta/phi*exp(-theta*t)*sin(phi*t);
        A2 = exp(-theta*t)*sin(phi*t)*1/phi;  
        Bb = 2*tau*t-2/kappa*(1-kappa*tau^2)*(exp(-theta*t)*cos(phi*t)-1)...
         -1/phi*tau*(3-kappa*tau^2)*exp(-theta*t)*sin(phi*t);
        Bc = -tau*(exp(-theta*t)*cos(phi*t)-1)- theta*tau/phi*exp(-theta*t)*sin(phi*t);
        B4= Bb*b(4)+Bc*c(4);

        Delta4_inf =A1*Delta3_inf + A2*Delta3_dot_inf+B4;
        Deltas(i_vec) =  Delta4_inf ;
        
    end
    figure(33+iter);clf;hold on;
    grid on
    set(gcf,'Position', [100+300*(iter-1),380,250,200]);
    %set(gca,'InnerPosition',[0.1 0.2 .9 .9]);
    plot(rs,Deltas./(a_bar),'linewidth',2)
    xlabel('$\bar{a}/\bar{v}$','interpreter','latex')
    ylabel('$\Delta/\bar{v}$','interpreter','latex')
    %zlabel('$\Delta$','interpreter','latex')
    %zlim([0,30]);
    title([subfig_ids{iter},' $\tau=$',num2str(tau),', $\kappa=$',num2str(kappa)],'interpreter','latex');
    saveas(gcf,['../../Figures/Delta_a',num2str(iter)],'png');
    
   
end

