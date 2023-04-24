clear;
digitsOld = digits(64);
kappa =vpa(3);
lambda = vpa(0.000001);
tau =vpa(.5);

ws = [0:0.05:3];
TFs_w = (kappa+i*lambda.*ws)./(-ws.^2+kappa+i*(lambda+kappa*tau).*ws);
TFs = abs(TFs_w);

omega = ws;
if lambda==0
   TF2s = kappa^2./(omega.^4+kappa*(kappa*tau^2-2).*omega.^2+kappa^2);
else
    
    A1 = lambda^2.*omega.^2+kappa^2;
    A2 = -2*kappa^2-kappa*(2-kappa*tau^2)*lambda^2+lambda^4+2*kappa*lambda^(3)*tau;
    A3 =  kappa^3*(kappa+kappa*lambda*tau+2*lambda^2)*(1-tau*lambda)*(lambda^2.*omega.^2+kappa^2).^(-1);
    Temp_As = A1+A2+A3;
    TF2s =lambda^4./(A1+A2+A3);  
    
end

check1 = kappa*tau^2+2*lambda*tau
if tau*lambda>=1
    w_star = 0;
    TF_star = 1;
elseif check1>=2 
    w_star = 0;
    TF_star = 1;
else
   
    if lambda==0
        w_star= sqrt(kappa*(2-kappa*tau^2)/2);
        PDF_star = 4/(4*kappa*tau^2-kappa^2*tau^4);
    else  
          Temp =  kappa^4+kappa^3*lambda^2*(2-kappa*tau^2-2*lambda*tau);
        w_star= sqrt((sqrt(Temp)-kappa^2))*lambda^(-1);
        
        %D = A2;
        %PDF_star = (2*lambda^(-4)*sqrt(Temp)+D)^-1;
        %D1= - lambda^(-4)*kappa^-2*Temp+1-kappa^2*lambda^-4;
        %
        %PDF_diff =  PDF_star - PDF_star1
        Temp1 =  Temp/kappa^2*lambda^-4;
        PDF_star = (1-(sqrt(Temp1)-kappa*lambda^-2)^2)^-1;

        TF_star = sqrt(PDF_star);
    end
   
    
    TF_star = sqrt(PDF_star);
    
end


figure(1);clf; hold all;
plot(ws, TFs,'b','linewidth',2);
plot(ws, sqrt(TF2s),':','linewidth',2);
plot(w_star,TF_star,'o','linewidth',4);
w_star
 TF_star

 
 digits(digitsOld)
