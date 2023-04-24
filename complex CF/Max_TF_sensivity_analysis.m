clear;
kappas =[0:0.05:3];
lambdas = [0:0.02:3];
[kappas_grid,lambdas_grid] = meshgrid(kappas,lambdas);
tau = .5;

w_stars = zeros(size(kappas_grid));
TF_stars = zeros(size(kappas_grid));
D_stars = zeros(size(kappas_grid));

kappa_stars = zeros(size(lambdas));

for i = 1:size(kappas_grid,1)
    kappa_stars(i)=nan;
    max_TF = -inf;
    for j = 1: size(kappas_grid,2)
        kappa = kappas_grid(i,j);
        lambda = lambdas_grid(i,j);
        check1 = kappa*tau^2+2*lambda*tau;
        if tau*lambda>=1
            w_star = 0;
            TF_star = 1;
            D=nan;
        elseif check1>=2 || lambda*kappa==0
            w_star = 0;
            TF_star = 1;
            D=nan;
        else
            D = -2.*kappa.^2.*lambda.^(-4)-2.*kappa.*lambda.^(-2)+1+2.*kappa.*lambda.^(-1).*tau+kappa.^2.*lambda.^(-2).*tau.^2;
            Temp =  kappa.^4+kappa.^3*lambda.^2.*(2-kappa.*tau.^2-2.*lambda.*tau);
            w_star= sqrt((sqrt(Temp)-kappa.^2)).*lambda.^(-1);
            PDF_star = (2.*lambda.^(-4)*sqrt(Temp)+D).^-1;
            TF_star = sqrt(PDF_star);
        end
        w_stars(i,j) =w_star;
        TF_stars(i,j) = TF_star;
        D_stars(i,j)=D;
        if TF_star > max_TF
            max_TF = TF_star;
            kappa_stars(i) = kappa;
        end
       
    end 
end

kappa_stars_sol = zeros(size(lambdas));
for i = 1:length(lambdas)
    lambda = lambdas(i);
    kappa_stars_sol(i) = max((1-lambda*tau)/(tau^2+lambda^-1*tau),0);
end

figure(1);clf; hold all;
surf(kappas,lambdas,w_stars)
%shading interp
xlabel('kappa')
ylabel('lambda')
zlabel('w^*')
view(3)
figure(2);clf; hold all;
surf(kappas,lambdas,TF_stars)

%shading interp
xlabel('kappa')
ylabel('lambda')
zlabel('TF^*')
view(3)

figure(3);clf; hold all;
plot(lambdas,kappa_stars_sol,'g')
plot(lambdas,kappa_stars,'g')


