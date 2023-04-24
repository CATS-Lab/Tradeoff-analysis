k = [0.0997,0.0693]';
t = [0.665,2.093]';
s0 = [13.558,7.730]';
l = [0.0926,0.0432]';
w_star = sqrt(4*k-2*k.^2.*t.^2)/2 ;
TF_star = sqrt(4./(4*k.*t.^2-k.^2.*t.^4));

w = [0:0.001:pi];
%T = (k+1i*l.*w)./(-w.^2+k+1i*(l+k.*t).*w);
T = (k)./(-w.^2+k+1i*k.*t.*w);

figure(1);cla;hold all;
plot(w,abs(T(1,:)),'r');
plot([w_star(1),w_star(1)],[0,5],'k:');
plot([w(1),w(end)],[TF_star(1),TF_star(1)],'k:');
plot(w,abs(T(2,:)),'b');

plot([w_star(2),w_star(2)],[0,5],'k:');
plot([w(1),w(end)],[TF_star(2),TF_star(2)],'k:');
xlabel('\omega (Hz)')
ylabel('T')

figure(2);cla;hold all;
plot(2*pi./w,abs(T(1,:)),'r');

plot(2*pi./w,abs(T(2,:)),'b');
xlim([0,50])
xlabel('period (sec)')
ylabel('T')
