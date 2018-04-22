% tst_ifsc - CHAPRO demonstration of GHA processing
function tst_ifsc
x=audioread('test/cat.wav')/5;
load('test/tst_ifsc')
y=wave;
figure(1); clf
nx=length(x);
ny=length(y);
tx=linspace(0,(nx - 1) / rate, nx);
ty=linspace(0,(ny - 1) / rate, ny);
plot(ty,y,'r',tx,x,'b')
xlabel('time (s)')
tlim=[min(ty) max(ty)]*1.05;
xylim=[min(min(x),min(y)) max(max(x),max(y))]*1.05;
axis([tlim xylim])
legend('output','input')
title('CHAPRO demonstration of GHA processing')
return
