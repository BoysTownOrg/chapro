% tst_ffsc - CHAPRO demonstration of GHA processing
function tst_ffsc
x=audioread('test/cat.wav')/5;
load('test/tst_ffsc')
y=wave;
figure(1); clf
n=length(wave);
t=linspace(0,(n - 1) / rate, n);
plot(t,y,'r',t,x,'b')
xlabel('time (s)')
ylim([-0.09 0.09])
legend('output','input')
title('CHAPRO demonstration of GHA processing')

return

