
rs = [1:1:20];
kappas = [0.1:0.1:2];

[r_mat,k_mat] = meshgrid(rs,kappas);
tau_underline_mat = zeros(size(k_mat));
g_underline_mat = zeros(size(k_mat));
Delta_v_underline_mat = zeros(size(k_mat));
 
a_bar = 1;

for i_mat = 1:size(k_mat,1)
    for j_mat = 1:size(k_mat,2)

        v_bar = r_mat(i_mat,j_mat)*a_bar;
        kappa = k_mat(i_mat,j_mat);
        tau = sqrt(2/kappa);
        
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
            Delta=  Delta4_inf ;

      
        g_underline_mat(i_mat,j_mat) = Delta/v_bar+tau;
        tau_underline_mat(i_mat,j_mat) = tau;
        Delta_v_underline_mat(i_mat,j_mat) =Delta/v_bar;
    end
end


ang1 =45;
ang2 = 30;


f1 = figure(60);clf;hold on;
grid on;
set(gcf,'Position', [100,380,250,200]);
surf(r_mat,k_mat,g_underline_mat);
 view(ang1,ang2)
ylabel('$\kappa$','interpreter','latex')
xlabel('$\bar{v}/\bar{a}$','interpreter','latex')
zlabel('$\underline{g}$','interpreter','latex')
title(['(a)'],'interpreter','latex');

f2 = figure(61);clf;hold on;
grid on;
set(gcf,'Position', [400,380,250,200]);
surf(r_mat,k_mat,g_underline_mat-g_star_mat);
 view(ang1,ang2)
ylabel('$\kappa$','interpreter','latex')
xlabel('$\bar{v}/\bar{a}$','interpreter','latex')
zlabel('$\underline{g}-g^*$','interpreter','latex')
title(['(b)'],'interpreter','latex');


f3 = figure(62);clf;hold on;
grid on;
set(gcf,'Position', [700,380,250,200]);
surf(r_mat,k_mat,Delta_v_underline_mat);
  view(ang1,ang2)
ylabel('$\kappa$','interpreter','latex')
xlabel('$\bar{v}/\bar{a}$','interpreter','latex')
zlabel('$\underline{\Delta}/\bar{v}$','interpreter','latex')
title(['(e)'],'interpreter','latex');

f4 = figure(63);clf;hold on;
grid on;
set(gcf,'Position', [1000,380,250,200]);
surf(r_mat,k_mat,Delta_v_underline_mat-Delta_v_star_mat);
 view(ang1,ang2)
ylabel('$\kappa$','interpreter','latex')
xlabel('$\bar{v}/\bar{a}$','interpreter','latex')
zlabel('$\underline{\Delta}/\bar{v}-\Delta^*/\bar{v}$','interpreter','latex')
title(['(f)'],'interpreter','latex');

f5 = figure(64);clf;hold on;
grid on;
set(gcf,'Position', [100,100,250,200]);
surf(r_mat,k_mat,tau_underline_mat);
 view(ang1,ang2)
ylabel('$\kappa$','interpreter','latex')
xlabel('$\bar{v}/\bar{a}$','interpreter','latex')
zlabel('$\underline{\tau}$','interpreter','latex')
title(['(c)'],'interpreter','latex');


f6 = figure(65);clf;hold on;
grid on;
set(gcf,'Position', [400,100,250,200]);
surf(r_mat,k_mat,tau_underline_mat-tau_star_mat);
 view(ang1,ang2)
ylabel('$\kappa$','interpreter','latex')
xlabel('$\bar{v}/\bar{a}$','interpreter','latex')
zlabel('$\underline{\tau}-\tau^*$','interpreter','latex')
title(['(d)'],'interpreter','latex');


% 
% saveas(f1,['../../Figures/g_underline'],'png');
% saveas(f2,['../../Figures/g_underline_minus_star'],'png');
% saveas(f3,['../../Figures/Delta_v_underline'],'png');
% saveas(f4,['../../Figures/Delta_v_underline_minus_star'],'png');
% saveas(f5,['../../Figures/tau_underline'],'png');
% saveas(f6,['../../Figures/tau_underline_minus_star'],'png');



