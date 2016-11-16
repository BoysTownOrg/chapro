% test FIR filterbank analysis
function tst_ffa
load('test/ffa_impulse')
y=real(y);
t=1000*(0:(length(y)-1))/rate;
figure(1);clf
plot(t,y)
title('FIR filterbank impulse responses')
y_lim=[min(min(y)) max(max(y))]* 1.05;
t_lim=[-0.5 15.5];
axis([t_lim y_lim])
xlabel('time (ms)')
figure(2);clf
H=ffa(y);
f=linspace(0,rate/2000,length(H))';
d=5*ones(size(f));
fm=max(max(f));
m_lim=[0.05 fm -50 10];
d_lim=[0.05 fm   0  8];
M=db(H);
D=gd(H,f);
D(M<-40)=NaN;
subplot(2,1,1)
semilogx(f,M)
axis(m_lim)
ylabel('magnitude (dB)')
title('FIR filterbank transfer functions')
subplot(2,1,2)
semilogx(f,D,f,d,':k')
axis(d_lim)
xlabel('frequency (kHz)')
ylabel('delay (ms)')
return

function y=db(x)
y=20*log10(max(eps,abs(x)));
return

function d=gd(H,f)
[nf,nc]=size(H);
d=zeros(nf,nc);
for k=1:nc
   p=unwrap(angle(H(:,k)))/(2*pi);
   d(:,k)=-cdif(p)./cdif(f);
end
return

function dx=cdif(x)
n=length(x);
dx=zeros(size(x));
dx(1)=x(2)-x(1);
dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2;
dx(n)=x(n)-x(n-1);
return
