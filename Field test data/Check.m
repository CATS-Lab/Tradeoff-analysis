load('data-Li')
i0=1;
i1=500;
v = zeros(i1,4);
 v(:,1)= highspeed_headway1(i0:i1,4);
 v(:,2)= highspeed_headway1(i0:i1,7);
 v(:,3)= highspeed_headway4(i0:i1,4);
 v(:,4)= highspeed_headway4(i0:i1,7);

figure(1);cla; hold all;
plot( v(:,1),'b:');
%plot( v(:,2),'b');
plot(v(:,3),'r:');
%plot(v(:,4),'r');


for i =1:4
    v(:,i)= v(:,i) - mean(v(:,i));
end
figure(2);

subplot(2,1,1);cla; hold all;
V=abs(fft(v));
is = [1:20];
xs = (i1-i0)./is;
plot(xs, V(is,1),'k');
plot(xs, V(is,2),'r');
subplot(2,1,2);cla; hold all;
plot(xs,V(is,3),'k');
plot(xs,V(is,4),'r');
ss = std(v);
A1 = ss(2)/ss(1)
A2 = ss(4)/ss(3)
