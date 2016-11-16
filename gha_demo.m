% gha_demo - CHAPRO demonstration of GHA processing
function gha_demo
load('test/gha_demo')

figure(2); clf
plot(t,y,'r',t,x,'b')
xlabel('time (s)')
ylim([-0.09 0.09])
legend('output','input')
title('CHAPRO demonstration of GHA processing')

return

