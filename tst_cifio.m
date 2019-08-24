% test gammatone filterbank input/output
function tst_cifio
load('test/gfio_impulse')
p_lim=[-0.1 1.1];
t_lim=[-1 30];
t=1000*(0:(length(y)-1))/rate;
figure(1);clf
subplot(1,2,1)
plot(t,x)
title('input')
axis([t_lim p_lim])
xlabel('time (ms)')
subplot(1,2,2)
plot(t,y)
title('output')
axis([t_lim p_lim])
xlabel('time (ms)')
% transfer function
figure(2);clf
dc=2.96;       % desired group delay (ms)
H=ffa(y)./ffa(x);
f=linspace(0,rate/2000,length(H))';
fm=max(f);
m_lim=[0.05 fm -5 5];
d_lim=[0.05 fm 1 9];
M=db(H);
D=gd(H,f);
d=dc*ones(size(f));
subplot(2,1,1)
semilogx(f,M)
axis(m_lim)
xlabel('frequency (kHz)')
ylabel('magnitude (dB)')
subplot(2,1,2)
semilogx(f,d,':',f,D)
axis(d_lim)
xlabel('frequency (kHz)')
ylabel('delay (ms)')
load('test/gfio_tone')
p_lim=[-2 2];
figure(3);clf
subplot(2,1,1)
plot(t,x)
title('input')
axis([t_lim p_lim])
subplot(2,1,2)
plot(t,y)
title('output')
axis([t_lim p_lim])
xlabel('time (ms)')
return

function y=db(x)
y=20*log10(max(eps,abs(x)));
return

function d=gd(x,f)
p=unwrap(angle(x))/(2*pi);
d=-cdif(p)./cdif(f);
return

function dx=cdif(x)
n=length(x);
dx=zeros(size(x));
dx(1)=x(2)-x(1);
dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2;
dx(n)=x(n)-x(n-1);
return

